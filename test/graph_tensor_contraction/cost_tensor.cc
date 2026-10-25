#include "../common/device.h"
#include "tnbp/framework/graph_contraction.h"
#include "tnbp/framework/graph_contraction_cost.h"
#include <complex>
#include <iostream>
#include <numeric>
using Scalar=std::complex<double>;
#ifdef TNBP_TEST_CUDA
using Tensor=tcapi::cuda::Tensor<Scalar>;
#else
using Tensor=gqten::tensor<Scalar>;
#endif
using Context=tcapi::context_handle_t<Tensor>;
using Edges=std::vector<std::pair<int,int>>;
using Shapes=std::map<int,std::vector<std::uint64_t>>;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void compare(Context& ctx,const char* name,const Edges& original,const std::vector<std::uint64_t>& dims,
             const Edges& schedule,const std::map<int,int>& priority,const std::map<int,int>& labels,MPI_Comm comm){
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    Shapes initial;
    for(auto site:tnbp::GetSiteIndexFromBond(original))
        for(auto e:tnbp::GetSurroundingBondIndex(site,original))initial[site].push_back(dims.at(e));
    auto estimated_edges=original;auto estimated_shapes=initial;
    std::vector<std::uint64_t> costs,elements;std::vector<std::array<int,3>> steps;
    std::vector<Shapes> shapes_trace;std::vector<Edges> edges_trace;
    const auto total=tnbp::graph_tensor_network_contraction_cost(estimated_edges,estimated_shapes,schedule,priority,labels,
        &costs,&elements,&steps,&shapes_trace,&edges_trace);
    require(costs.size()==steps.size() && total==std::accumulate(costs.begin(),costs.end(),std::uint64_t(0)),"total MAC count");
    for(bool empty_root:{false,true}){
        Edges edges=original;std::map<int,int> owners;std::vector<int> ids;std::vector<Tensor> tensors;
        auto initialize=[&]{
            owners.clear();ids.clear();tensors.clear();
            for(const auto& item:initial){
                owners[item.first]=empty_root?size-1:(item.first+1)%size;
                if(owners[item.first]!=rank)continue;
                ids.push_back(item.first);tcapi::shape_t<Tensor> shape(item.second.begin(),item.second.end());
                tensors.push_back(tcapi::zeros<Tensor>(ctx,shape));
                tcapi::set_elem(ctx,tensors.back(),tcapi::elem_coors_t<Tensor>(shape.size(),0),Scalar(1));
            }
        };
        auto check=[&](const Shapes& shapes){
            require(ids.size()==tensors.size() && owners.size()==shapes.size(),"tensor/ID/owner count");
            std::size_t local_count=0;for(const auto& item:shapes)if(owners.at(item.first)==rank)++local_count;
            require(ids.size()==local_count,"unexpected local tensor");
            for(const auto& item:shapes){
                auto it=std::find(ids.begin(),ids.end(),item.first);
                require((it!=ids.end())==(owners.at(item.first)==rank),"owner-local result");
                if(it==ids.end())continue;
                const auto& tensor=tensors.at(it-ids.begin());auto actual=tcapi::shape(ctx,tensor);
                require(std::vector<std::uint64_t>(actual.begin(),actual.end())==item.second,"actual shape differs from cost shape (including external neighbors)");
                require(tcapi::get_elem(ctx,tensor,tcapi::elem_coors_t<Tensor>(actual.size(),0))==Scalar(1),"sparse fixture contraction value");
            }
        };
        initialize();
        for(std::size_t i=0;i<steps.size();++i){
            const auto op=steps[i];const auto expected_owner=owners.at(op[0]);
            tnbp::graph_tensor_contraction(ctx,edges,tensors,ids,owners,op[0],op[1],op[2],comm);
            require(edges==edges_trace[i] && owners.at(op[2])==expected_owner,"step edges or owner differs");
            check(shapes_trace[i]);
            if(owners.at(op[2])==rank){
                const auto index=std::find(ids.begin(),ids.end(),op[2])-ids.begin();std::uint64_t n=1;
                for(auto d:tcapi::shape(ctx,tensors[index]))n*=d;
                require(n==elements[i],"actual output size differs from cost size");
            }
        }
        require(edges==estimated_edges,"step final edges");
        // Exercise the real network executor independently of the returned steps.
        // This catches drift in priority/label selection or schedule rewriting.
        edges=original;initialize();
        tnbp::graph_tensor_network_contraction(ctx,edges,tensors,ids,owners,schedule,priority,labels,comm);
        require(edges==estimated_edges,"network executor final edges");check(estimated_shapes);
    }
    if(rank==0)std::cout<<name<<": actual tensor steps/shapes/owners and network schedule PASS; MACs="<<total<<'\n';
}
int main(int argc,char** argv){
    MPI_Init(&argc,&argv);int rank;MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    try{
        initialize_tensor_device(MPI_COMM_WORLD);Context ctx;tcapi::create_context(ctx);
        const Edges triangle{{0,1},{0,2},{1,2}};
        compare(ctx,"triangle",triangle,{2,3,5},triangle,{}, {},MPI_COMM_WORLD);
        compare(ctx,"reverse tie",triangle,{2,3,5},Edges{{2,1},{0,1},{0,2}},{}, {},MPI_COMM_WORLD);
        compare(ctx,"priority/fresh label",triangle,{2,3,5},triangle,{{1,5}},{{1,-9}},MPI_COMM_WORLD);
        Edges k4{{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
        compare(ctx,"K4 multiple neighbors",k4,{2,3,2,3,2,3},k4,{}, {},MPI_COMM_WORLD);
        std::reverse(k4.begin(),k4.end());for(auto& edge:k4)std::swap(edge.first,edge.second);
        compare(ctx,"K4 reversed edge order",k4,{3,2,3,2,3,2},k4,{{3,1}}, {},MPI_COMM_WORLD);
        Edges square{{0,1},{0,2},{0,3},{1,2},{2,3},{1,4},{2,5},{3,6}};
        compare(ctx,"square partial",square,{3,2,3,2,2,3,2,3},{{1,0},{2,1},{2,3}},{{1,1},{2,2}}, {},MPI_COMM_WORLD);
        Edges honeycomb;
        for(int row=0;row<2;++row)for(int cell=0;cell<3;++cell){int a=6*row+2*cell;honeycomb.emplace_back(a,a+1);if(cell)honeycomb.emplace_back(a,a-1);if(row)honeycomb.emplace_back(a,a-5);}
        compare(ctx,"honeycomb full",honeycomb,std::vector<std::uint64_t>(honeycomb.size(),2),honeycomb,{}, {},MPI_COMM_WORLD);
        compare(ctx,"empty schedule",square,{3,2,3,2,2,3,2,3},{},{}, {},MPI_COMM_WORLD);
        tcapi::destroy_context(ctx);if(rank==0)std::cout<<"graph cost versus numerical execution PASS\n";
    }catch(const std::exception& e){std::cerr<<"rank "<<rank<<": "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}
    MPI_Finalize();return 0;
}
