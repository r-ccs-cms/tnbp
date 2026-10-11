#include "tnbp/ptns/function/bmps.h"
#include "../common/dense_fixture.h"
#include <iostream>
using Tensor=gqten::tensor<Scalar>;
using Context=tcapi::context_handle_t<Tensor>;
using Edges=std::vector<std::pair<int,int>>;
// Compare the actual two-line tensor subset used by sampling, with and without
// contracting the far upper line. Diagonals exercise common-neighbor bond fusion.
void upper_independence(Context& ctx,MPI_Comm comm) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    for(bool reversed:{false,true}) for(bool manual:{false,true}) {
        Edges edges{{0,1},{2,3},{4,5},{6,7},
                    {0,2},{1,3},{0,3},{2,4},{3,5},{2,5},
                    {4,6},{5,7},{4,7}};
        if(reversed) std::reverse(edges.begin(),edges.end());
        std::vector<Tensor> all;std::vector<int> ids;std::map<int,int> owners;
        for(int site=0;site<8;++site) {
            owners[site]=site%size;
            if(owners[site]!=rank) continue;
            ids.push_back(site);
            const auto degree=tnbp::GetSurroundingBondIndex(site,edges).size();
            auto t=tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>(degree,2));
            for(int bits=0;bits<(1<<degree);++bits) {
                tcapi::elem_coors_t<Tensor> c;std::vector<int> values;
                for(std::size_t k=0;k<degree;++k) {
                    c.push_back((bits>>k)&1);values.push_back((bits>>k)&1);
                }
                tcapi::set_elem(ctx,t,c,entry(site,site%2,values));
            }
            all.push_back(std::move(t));
        }
        auto compact_edges=edges;auto compact_ids=ids;auto compact_owners=owners;
        std::vector<Tensor> compact;
        for(const auto& t:all) compact.push_back(tcapi::copy(ctx,t));
        const std::vector<std::vector<int>> upper_groups{{4,6},{5,7}};
        tnbp::bmps_increment_line(ctx,std::vector<int>{6,7},std::vector<int>{4,5},
            compact_edges,compact,compact_ids,compact_owners,comm,
            manual ? &upper_groups : nullptr);
        // Discard every upper tensor. Both lower updates must work with only 0..3.
        auto subset=[&](const std::vector<Tensor>& source,const std::vector<int>& labels) {
            std::vector<Tensor> result;
            for(std::size_t i=0;i<labels.size();++i)
                if(labels[i]<4) result.push_back(tcapi::copy(ctx,source[i]));
            return result;
        };
        auto full= subset(all,ids), packed=subset(compact,compact_ids);
        std::vector<int> full_ids,packed_ids;
        for(int site:ids) if(site<4) full_ids.push_back(site);
        for(int site:compact_ids) if(site<4) packed_ids.push_back(site);
        const std::vector<int> input{0,1},output{2,3};
        const std::vector<std::vector<int>> groups{{2,0,4},{3,1,5}};
        const auto* plan=manual ? &groups : nullptr;
        if(tnbp::bmps_absorption_schedule(input,output,edges,plan)!=
           tnbp::bmps_absorption_schedule(input,output,compact_edges,plan))
            throw std::runtime_error("upper contraction changes lower schedule");
        tnbp::bmps_increment_line(ctx,input,output,edges,full,full_ids,owners,comm,plan);
        tnbp::bmps_increment_line(ctx,input,output,compact_edges,packed,packed_ids,compact_owners,comm,plan);
        auto compare=[&] {
            if(full_ids!=packed_ids) throw std::runtime_error("lower site labels differ");
            for(std::size_t i=0;i<full.size();++i) {
                const auto shape=tcapi::shape(ctx,full[i]);
                if(shape!=tcapi::shape(ctx,packed[i])) throw std::runtime_error("lower shapes differ");
                int count=1;for(auto d:shape) count*=d;
                for(int n=0;n<count;++n) {
                    int v=n;tcapi::elem_coors_t<Tensor> c;
                    for(auto d:shape) {c.push_back(v%d);v/=d;}
                    const auto a=tcapi::get_elem(ctx,full[i],c),b=tcapi::get_elem(ctx,packed[i],c);
                    if(std::abs(a-b)>1.e-11*(1.+std::abs(a)))
                        throw std::runtime_error("upper contraction changes lower tensor");
                }
            }
            for(int site:output) {
                if(owners.at(site)!=site%size || compact_owners.at(site)!=site%size)
                    throw std::runtime_error("lower owner changed");
                Edges a,b;
                for(auto e:edges) if(e.first==site||e.second==site) a.push_back(e);
                for(auto e:compact_edges) if(e.first==site||e.second==site) b.push_back(e);
                if(a!=b) throw std::runtime_error("lower incident edge order changed");
            }
        };
        compare();
        // Force truncation (bond dimension 1), and compare after BP as well.
        auto truncate=[&](const Edges& graph,std::vector<Tensor>& ts,
                          const std::vector<int>& labels,const std::map<int,int>& ranks) {
            std::vector<Tensor> messages;std::vector<int> message_ids;
            std::vector<tcapi::bond_dim_t<Tensor>> dims;std::vector<double> errors;
            tnbp::bmps_bp_truncation(ctx,graph,output,ts,labels,ranks,messages,message_ids,
                30,1.e-12,tcapi::bond_dim_t<Tensor>(1),1.e-14,1.e-14,dims,errors,comm);
        };
        truncate(edges,full,full_ids,owners);
        truncate(compact_edges,packed,packed_ids,compact_owners);
        compare();
    }
    if(rank==0) std::cout<<"upper contraction independence: auto/manual, forward/reversed edges, absorption/BP truncation PASS"<<std::endl;
}
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    try {
        Context ctx;tcapi::create_context(ctx);
        upper_independence(ctx,MPI_COMM_WORLD);
        tcapi::destroy_context(ctx);
    } catch(const std::exception& e) {std::cerr<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,1);}
    MPI_Finalize();
}
