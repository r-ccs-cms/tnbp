#pragma once
#include "tnbp/framework/graph_contraction.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <string>

using Edges=std::vector<std::pair<int,int>>;
#ifdef TNBP_TEST_CUDA
template<class S> using TensorOf=tcapi::cuda::Tensor<S>;
#else
template<class S> using TensorOf=gqten::tensor<S>;
#endif

template<class S> S value(std::size_t site,std::size_t index) {
    const double x=0.2+std::sin(0.71*(index+1)*(site+1));
    if constexpr(std::is_same_v<S,double>) return x;
    else return S(x,0.13*std::cos(0.37*(index+2)*(site+1)));
}

// Dense enumeration of internal bonds; leave every external bond open. This
// reference does not use TNBP graph contraction or the tested factorization.
template<class T>
std::vector<tcapi::elem_t<T>> boundary(tcapi::context_handle_t<T>& ctx,
    const std::vector<int>& line,const Edges& edges,const std::vector<T>& local_t,
    const std::vector<int>& ids,const std::map<int,int>& owners,MPI_Comm comm) {
    using S=tcapi::elem_t<T>;
    int rank;MPI_Comm_rank(comm,&rank);
    std::set<int> members(line.begin(),line.end());
    std::vector<std::vector<S>> data;
    std::vector<std::vector<std::size_t>> bonds;
    std::vector<std::size_t> dims(edges.size(),0);
    for(auto site:line) {
        T t;
        if(rank==owners.at(site)) {
            auto i=std::find(ids.begin(),ids.end(),site)-ids.begin();
            t=tcapi::copy(ctx,local_t.at(i));
        }
        tnbp::MpiBcast(ctx,t,owners.at(site),comm);
        const auto shape=tcapi::shape(ctx,t);
        const auto incident=tnbp::GetSurroundingBondIndex(site,edges);
        if(shape.size()!=incident.size()) throw std::runtime_error("reference rank mismatch");
        bonds.emplace_back(incident.begin(),incident.end());
        for(std::size_t k=0;k<shape.size();++k) {
            auto& d=dims[incident[k]];
            if(d && d!=static_cast<std::size_t>(shape[k])) throw std::runtime_error("bond dimensions mismatch");
            d=shape[k];
        }
        std::vector<S> flat(tcapi::size(ctx,t));
        tcapi::to_range(ctx,t,flat.begin(),[shape](const auto& c){return tnbp::address_from_coor(shape,c);});
        data.push_back(std::move(flat));
    }
    std::vector<std::size_t> internal,external;
    for(std::size_t e=0;e<edges.size();++e) {
        const int count=members.count(edges[e].first)+members.count(edges[e].second);
        if(count==2) internal.push_back(e);
        else if(count==1) external.push_back(e);
    }
    // External edge order must remain unchanged, so this is a common output basis.
    std::size_t ni=1,no=1;
    for(auto e:internal) ni*=dims[e];
    for(auto e:external) no*=dims[e];
    if(ni*no>5000000) throw std::runtime_error("dense reference too large");
    std::vector<S> out(no,S{});
    std::vector<std::size_t> coords(edges.size());
    for(std::size_t o=0;o<no;++o) {
        auto v=o;for(auto e:external) {coords[e]=v%dims[e];v/=dims[e];}
        for(std::size_t i=0;i<ni;++i) {
            v=i;for(auto e:internal) {coords[e]=v%dims[e];v/=dims[e];}
            S product=1;
            for(std::size_t s=0;s<line.size();++s) {
                std::size_t addr=0,stride=1;
                for(auto e:bonds[s]) {addr+=coords[e]*stride;stride*=dims[e];}
                product*=data[s][addr];
            }
            out[o]+=product;
        }
    }
    return out;
}

template<class S, class ShiftEdges, class ShiftTensors>
void run_shift_fixture(tcapi::context_handle_t<TensorOf<S>>& ctx,
         const std::string& name,Edges edges,const std::vector<int>& line,
         int ownership_shift,bool zero,MPI_Comm comm,ShiftEdges shift_edges,ShiftTensors shift_tensors,bool single=false) {
    using T=TensorOf<S>;
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    std::map<int,int> owners;
    std::vector<int> ids;
    std::vector<T> ts;
    // A scalar spectator has no edge and must not be affected by temporary labels.
    const int spectator=std::numeric_limits<int>::max();owners[spectator]=0;
    if(rank==0) {ids.push_back(spectator);ts.push_back(tcapi::zeros<T>(ctx,{}));tcapi::set_elem(ctx,ts.back(),{},S(7));}
    for(std::size_t s=0;s<line.size();++s) {
        const int site=line[s],owner=(s+ownership_shift)%size;owners[site]=owner;
        if(owner!=rank) continue;
        ids.push_back(site);
        const auto bonds=tnbp::GetSurroundingBondIndex(site,edges);
        tcapi::shape_t<T> shape;
        for(auto e:bonds) shape.push_back((e%3==0) ? 3 : 2);
        std::size_t count=1;for(auto d:shape) count*=d;
        std::vector<S> flat(count);
        for(std::size_t i=0;i<count;++i) flat[i]=zero ? S{} : value<S>(s,i);
        ts.push_back(tcapi::assign_from_range<T>(ctx,shape,flat.begin(),
            [shape](const auto& c){return tnbp::address_from_coor(shape,c);}));
    }
    const auto initial_ids=ids;const auto initial_owners=owners;
    const auto expected=boundary(ctx,line,edges,ts,ids,owners,comm);
    const auto old_edges=edges;
    auto planned_edges=edges;
    if(single) tnbp::graph_shift_bonds_to_neighbor_for_edges(planned_edges,line.back(),line[line.size()-2],std::vector<int>{line.front()});
    else shift_edges(line,planned_edges);
    if(single) {
        const int source=line.back(),destination=line[line.size()-2],other=line.front();
        std::set<int> required{source,destination};
        if(tnbp::find_edge(edges,destination,other)!=edges.end()) required.insert(other);
        std::vector<T> partial;std::vector<int> partial_ids;std::map<int,int> partial_owners;
        for(int site:required) partial_owners[site]=owners.at(site);
        for(std::size_t i=0;i<ids.size();++i) if(required.count(ids[i])) {
            partial_ids.push_back(ids[i]);partial.push_back(tcapi::copy(ctx,ts[i]));
        }
        const auto old_ids=partial_ids;const auto old_owners=partial_owners;
        tnbp::graph_shift_bond_to_neighbor(ctx,edges,partial,partial_ids,partial_owners,
                                         source,destination,other,comm);
        if(partial_ids!=old_ids || partial_owners!=old_owners)
            throw std::runtime_error("single shift changed metadata");
        for(std::size_t i=0;i<partial_ids.size();++i) {
            auto j=std::find(ids.begin(),ids.end(),partial_ids[i])-ids.begin();
            ts[j]=tcapi::copy(ctx,partial[i]);
        }
        std::set<std::pair<int,int>> want,got;
        for(auto e:old_edges) {
            if(tnbp::same_edge(e,tnbp::make_edge(source,other))) e=tnbp::make_edge(destination,other);
            want.insert(tnbp::make_edge(e.first,e.second));
        }
        want.insert(tnbp::make_edge(source,destination));
        for(auto e:edges) got.insert(tnbp::make_edge(e.first,e.second));
        if(want!=got) throw std::runtime_error("single shift changed unselected edges");
    } else shift_tensors(ctx,line,edges,ts,ids,owners,comm);
    if(ids!=initial_ids || owners!=initial_owners) throw std::runtime_error("site/owner metadata changed");
    if(rank==0 && tcapi::get_elem(ctx,ts[0],{})!=S(7)) throw std::runtime_error("spectator changed");
    std::map<int,std::size_t> pos;for(std::size_t i=0;i<line.size();++i) pos[line[i]]=i;
    if(edges!=planned_edges) throw std::runtime_error("ordered tensor/graph-only shift edges disagree");
    Edges old_external,new_external;
    for(auto e:old_edges) if(!pos.count(e.first)||!pos.count(e.second)) old_external.push_back(tnbp::make_edge(e.first,e.second));
    std::set<std::pair<int,int>> unique;
    for(auto e:edges) {
        if(!unique.insert(tnbp::make_edge(e.first,e.second)).second || e.first==e.second)
            throw std::runtime_error("invalid resulting graph");
        if(pos.count(e.first)&&pos.count(e.second)) {
            auto a=pos.at(e.first),b=pos.at(e.second);
            if(!single && (a>b ? a-b : b-a)!=1) throw std::runtime_error("long bond remains");
        } else new_external.push_back(tnbp::make_edge(e.first,e.second));
    }
    if(old_external!=new_external) throw std::runtime_error("external edges changed");
    const auto actual=boundary(ctx,line,edges,ts,ids,owners,comm);
    if(actual.size()!=expected.size()) throw std::runtime_error("boundary shape changed");
    double max_error=0,scale=0;
    for(std::size_t i=0;i<expected.size();++i) {
        if(!std::isfinite(std::abs(actual[i]))) throw std::runtime_error("nonfinite result");
        max_error=std::max(max_error,double(std::abs(actual[i]-expected[i])));
        scale=std::max(scale,double(std::abs(expected[i])));
    }
    if(max_error>1.e-10*(1+scale)) throw std::runtime_error(name+": boundary tensor changed");
    if(!single) {
    // A second pass must do nothing, including leaving tensor gauges unchanged.
    std::vector<std::vector<S>> before;
    for(auto& t:ts) {
        auto shape=tcapi::shape(ctx,t);std::vector<S> data(tcapi::size(ctx,t));
        tcapi::to_range(ctx,t,data.begin(),[shape](const auto& c){return tnbp::address_from_coor(shape,c);});
        before.push_back(std::move(data));
    }
    const auto chain_edges=edges;
    shift_tensors(ctx,line,edges,ts,ids,owners,comm);
    if(edges!=chain_edges) throw std::runtime_error("second pass changed graph");
    for(std::size_t k=0;k<ts.size();++k) {
        auto shape=tcapi::shape(ctx,ts[k]);std::vector<S> data(tcapi::size(ctx,ts[k]));
        tcapi::to_range(ctx,ts[k],data.begin(),[shape](const auto& c){return tnbp::address_from_coor(shape,c);});
        if(data!=before[k]) throw std::runtime_error("second pass changed tensor");
    }
    }
    if(rank==0) std::cout<<name<<" "<<(std::is_same_v<S,double> ? "real" : "complex")
                        <<" error="<<max_error<<" PASS"<<std::endl;
}
