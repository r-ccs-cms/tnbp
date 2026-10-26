#include "../common/device.h"
#include "tnbp/ptns/function/mcmc_sampling.h"
#include <array>
#include <iostream>
#include <random>
#include <stdexcept>

using Scalar=std::complex<double>;
#ifdef TNBP_TEST_CUDA
using Tensor=tcapi::cuda::Tensor<Scalar>;
#else
using Tensor=gqten::tensor<Scalar>;
#endif
using States=std::vector<std::vector<int>>;
using Samples=std::vector<States>;
using Lines=std::vector<std::vector<int>>;
using Context=tcapi::context_handle_t<Tensor>;
void require(bool value,const char* message) {if(!value) throw std::runtime_error(message);}

void run(Context& ctx,int rank,int size) {
    // Exactly one particle on a square. A cumulative charge on the path 0-1-3-2
    // gives amplitude weight[occupied site]; the remaining edge has dimension 1.
    const std::vector<std::pair<int,int>> edges{{0,1},{0,2},{1,3},{2,3}};
    const std::array<int,4> path{0,1,3,2};
    const std::array<Scalar,4> weight{Scalar(1,.2),Scalar(.8,-.3),Scalar(1.2,.1),Scalar(.9,.4)};
    std::vector<int> ids;
    std::map<int,int> owners;
    std::vector<Tensor> W;
    for(int site=0;site<4;++site) {
        owners[site]=site%std::min(size,2); // third rank deliberately has no tensors
        if(owners[site]!=rank) continue;
        ids.push_back(site);
        const auto pos=std::find(path.begin(),path.end(),site)-path.begin();
        tcapi::shape_t<Tensor> shape;
        std::vector<int> neighbors;
        for(const auto& edge:edges) if(edge.first==site || edge.second==site) {
            neighbors.push_back(edge.first==site ? edge.second : edge.first);
            shape.push_back(edge==std::make_pair(0,2) ? 1 : 2);
        }
        shape.push_back(2);
        auto t=tcapi::zeros<Tensor>(ctx,shape);
        for(int incoming=0;incoming<2;++incoming) for(int s=0;s<2;++s) {
            if(pos==0 && incoming!=0) continue;
            const int outgoing=incoming+s;
            if(outgoing>1 || (pos==3 && outgoing!=1)) continue;
            tcapi::elem_coors_t<Tensor> c(neighbors.size()+1,0);
            for(std::size_t a=0;a<neighbors.size();++a) {
                if(pos>0 && neighbors[a]==path[pos-1]) c[a]=incoming;
                if(pos<3 && neighbors[a]==path[pos+1]) c[a]=outgoing;
            }
            c.back()=s;
            tcapi::set_elem(ctx,t,c,s ? weight[site] : Scalar(1));
        }
        W.push_back(std::move(t));
    }
    const std::array<Lines,2> layouts{Lines{{0,1},{2,3}},Lines{{0,2},{1,3}}};
    // A manual / B automatic also exercises clearing layout-specific grouping.
    const std::array<tnbp::MCMCGrouping<int>,2> groups{
        tnbp::MCMCGrouping<int>{{{0,2},{1,3}},{{2,0},{3,1}}},tnbp::MCMCGrouping<int>{}};
    States initial(2,std::vector<int>(ids.size()));
    for(int b=0;b<2;++b) for(std::size_t i=0;i<ids.size();++i) initial[b][i]=(ids[i]==b*3);
    auto sample_phase=[&](const States& current,int phase,int count,std::mt19937& rng) {
        tnbp::MCMCWorkspace<Tensor,int> work;
        work.lines=layouts[phase];work.manual_grouping=groups[phase];
        tnbp::mcmc_initialize(ctx,edges,W,ids,owners,current,10,1e-12,16,1e-12,1e-12,work,MPI_COMM_WORLD);
        Samples out(count,States(current.size()));out[0]=current;
        tnbp::mcmc_sampling(ctx,edges,W,ids,owners,out,2,rng,10,1e-12,16,1e-12,1e-12,work,MPI_COMM_WORLD,2);
        return out;
    };
    auto combined=[&](const States& current,std::array<int,2> trips,Samples& out,std::mt19937& rng) {
        tnbp::mcmc_sampling_alternating(ctx,edges,W,ids,owners,current,layouts,groups,trips,
            out,rng,10,1e-12,16,1e-12,1e-12,MPI_COMM_WORLD);
    };
    for(const auto trips:{std::array<int,2>{1,1},std::array<int,2>{2,3}}) {
        std::mt19937 rng(123+rank),reference_rng=rng;
        Samples actual(4);combined(initial,trips,actual,rng);
        auto current=initial;
        for(const auto& output:actual) {
            for(int phase=0;phase<2;++phase) current=sample_phase(current,phase,trips[phase],reference_rng).back();
            require(output==current,"alternating state differs from explicit composition");
        }
        require(rng==reference_rng,"alternating RNG differs from explicit composition");
    }
    std::mt19937 rng(777+rank),saved=rng;
    Samples empty;combined(initial,{2,3},empty,rng);
    require(rng==saved && empty.empty(),"zero output changed RNG/output");
    auto global=[&](const std::vector<int>& local) {
        int here[4]={},all[4]={};
        for(std::size_t i=0;i<ids.size();++i) here[ids[i]]=local[i];
        MPI_Allreduce(here,all,4,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        return std::array<int,4>{all[0],all[1],all[2],all[3]};
    };
    auto horizontal=sample_phase(initial,0,32,rng);
    for(const auto& sample:horizontal) for(int b=0;b<2;++b) {
        auto s=global(sample[b]);
        require(s[0]+s[1]==1-b && s[2]+s[3]==b,"horizontal changed row charge");
    }
    Samples out(4096+128);combined(States{initial[0]},{1,1},out,rng);
    std::array<int,4> histogram{};
    for(std::size_t j=0;j<out.size();++j) {
        auto state=global(out[j][0]);int total=0;
        for(int s:state) total+=s;
        require(total==1,"alternating violated total charge");
        if(j>=128) for(int i=0;i<4;++i) histogram[i]+=state[i];
    }
    double norm=0,max_error=0,tv=0;
    for(auto w:weight) norm+=std::norm(w);
    for(int i=0;i<4;++i) {
        require(histogram[i]>0,"alternating failed to reach every site");
        const double error=std::abs(histogram[i]/4096.0-std::norm(weight[i])/norm);
        max_error=std::max(max_error,error);tv+=error/2;
    }
    require(max_error<.04 && tv<.10,"alternating Born distribution mismatch");
    if(rank==0) std::cout<<"alternating PASS: composition {1,1}/{2,3}, 2 batches, RNG, zero output, row/total charge; 4096 samples max_error="<<max_error<<" TV="<<tv<<'\n';
}
// Different line counts force a fresh three-line boundary in the second phase.
void rectangle(Context& ctx,int rank,int size) {
    const std::vector<std::pair<int,int>> edges{{0,1},{1,2},{3,4},{4,5},{0,3},{1,4},{2,5}};
    const std::array<Lines,2> layouts{Lines{{0,1,2},{3,4,5}},Lines{{0,3},{1,4},{2,5}}};
    std::array<tnbp::MCMCGrouping<int>,2> groups;
    groups[1]=tnbp::make_mcmc_grouping(edges,layouts[1]);
    std::map<int,int> owners;
    std::vector<int> ids;
    std::vector<Tensor> W;
    std::mt19937 source(42+rank);
    for(int site=0;site<6;++site) {
        owners[site]=site%size;
        if(owners[site]!=rank) continue;
        ids.push_back(site);
        tcapi::shape_t<Tensor> shape(tnbp::GetSurroundingBondIndex(site,edges).size()+1,2);
        W.push_back(tcapi::random<Tensor>(ctx,shape,source));
    }
    States initial(2,std::vector<int>(ids.size()));
    for(int b=0;b<2;++b) for(std::size_t i=0;i<ids.size();++i) initial[b][i]=(ids[i]+b)%2;
    std::mt19937 rng(789+rank),reference_rng=rng;
    Samples actual(3);
    tnbp::mcmc_sampling_alternating(ctx,edges,W,ids,owners,initial,layouts,groups,{2,3},
        actual,rng,10,1e-12,32,1e-12,1e-12,MPI_COMM_WORLD);
    auto current=initial;
    for(const auto& output:actual) {
        for(int phase=0;phase<2;++phase) {
            tnbp::MCMCWorkspace<Tensor,int> work;
            work.lines=layouts[phase];work.manual_grouping=groups[phase];
            tnbp::mcmc_initialize(ctx,edges,W,ids,owners,current,10,1e-12,32,1e-12,1e-12,work,MPI_COMM_WORLD);
            Samples step(phase==0 ? 2 : 3,States(current.size()));step[0]=current;
            tnbp::mcmc_sampling(ctx,edges,W,ids,owners,step,2,reference_rng,10,1e-12,32,1e-12,1e-12,work,MPI_COMM_WORLD,2);
            current=std::move(step.back());
        }
        require(output==current,"rectangle state differs from explicit composition");
    }
    require(rng==reference_rng,"rectangle RNG differs from explicit composition");
    if(rank==0) std::cout<<"rectangle PASS: 2/3 lines, automatic/manual groups, 2 batches, 3 cycles, {2,3} round trips\n";
}
void fatal(Context& ctx,int which) {
    // Empty local tensor arrays are intentional: invalid second-layout metadata
    // must be rejected before any first-phase numerical initialization.
    const std::vector<std::pair<int,int>> edges{{0,1},{0,2},{1,3},{2,3}};
    std::array<Lines,2> layouts{Lines{{0,1},{2,3}},Lines{{0,2},{1,3}}};
    std::array<tnbp::MCMCGrouping<int>,2> groups;
    std::array<int,2> trips{1,1};
    if(which==0) trips[1]=0;
    if(which==1) layouts[1]={{0,2},{1,1}};
    if(which==2) groups[1]={{{0},{2}},{{1},{3}}};
    std::mt19937 rng(1);Samples output(which==3 ? 0 : 1);
    if(which==3) layouts[1].clear();
    tnbp::mcmc_sampling_alternating(ctx,edges,std::vector<Tensor>{},std::vector<int>{},
        std::map<int,int>{{0,0},{1,0},{2,0},{3,0}},States(1),layouts,groups,trips,
        output,rng,10,1e-12,16,1e-12,1e-12,MPI_COMM_WORLD);
    throw std::runtime_error("invalid alternating request accepted");
}
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);initialize_tensor_device(MPI_COMM_WORLD);
    int rank,size;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
    try {
        Context ctx;tcapi::create_context(ctx);
        if(argc==3 && std::string(argv[1])=="fatal") fatal(ctx,std::stoi(argv[2]));
        else {run(ctx,rank,size);rectangle(ctx,rank,size);}
        tcapi::destroy_context(ctx);
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}
    MPI_Finalize();
}
