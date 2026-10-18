#include "tnbp/tnbp.h"
#include <complex>
#include <limits>
#include <iostream>
#ifdef TNBP_TEST_CUDA
#include <cuda_runtime.h>
template<class E> using Tensor=tcapi::cuda::Tensor<E>;
#else
template<class E> using Tensor=gqten::tensor<E>;
#endif
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
template<class E> void run(MPI_Comm comm){
  using T=Tensor<E>;using R=tnbp::real_t<T>;
  tnbp::context_handle_t<T> ctx;tcapi::create_context(ctx);
  int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
  const std::vector<std::pair<int,int>> edges{{0,1},{1,2}};
  const std::map<int,int> owners{{0,0},{1,1%size},{2,2%size}};
  std::vector<int> ids;
  for(int j=0;j<2;++j)if(owners.at(edges[j].first)==rank || owners.at(edges[j].second)==rank)ids.push_back(j);
  std::vector<T> old(2*ids.size()),next(2*ids.size());
  E phase;
  if constexpr(std::is_arithmetic_v<E>)phase=E(-1);
  else phase=std::polar(R(1),R(.7));
  for(std::size_t j=0;j<ids.size();++j)for(int direction=0;direction<2;++direction){
    const auto edge=edges[ids[j]];bool own=owners.at(direction?edge.second:edge.first)==rank;
    const auto slot=j+direction*ids.size();
    old[slot]=tcapi::zeros<T>(ctx,{2,2});next[slot]=tcapi::zeros<T>(ctx,{2,2});
    if(own){
      tcapi::for_each_with_coors(ctx,old[slot],[](auto& x,const auto& c){if(c[0]==0&&c[1]==0)x=E(3);});
      tcapi::for_each_with_coors(ctx,next[slot],[phase](auto& x,const auto& c){if(c[0]==1&&c[1]==1)x=E(5)*phase;});
    }
    // Non-owner replicas intentionally stay zero: they must never be normalized.
  }
  for(std::size_t j=0;j<ids.size();++j)for(int direction=0;direction<2;++direction){
    const auto edge=edges[ids[j]];
    if(owners.at(direction?edge.second:edge.first)!=rank)continue;
    const auto slot=j+direction*ids.size();
    tnbp::NormalizeMessage(ctx,old[slot],comm);
    tnbp::NormalizeMessage(ctx,next[slot],comm);
  }
  std::vector<E> before;
  for(auto* messages:{&old,&next})for(const auto& t:*messages)for(int i=0;i<2;++i)for(int j=0;j<2;++j)before.push_back(tcapi::get_elem(ctx,t,{i,j}));
  R residual=0;tnbp::BeliefPropagationDifference(ctx,edges,owners,old,next,ids,comm,residual);
  const R tol=R(128)*std::numeric_limits<R>::epsilon();
  require(std::abs(residual-std::sqrt(R(2)))<tol,"independent orthogonal-message distance / global directed-edge weighting");
  std::size_t at=0;
  for(auto* messages:{&old,&next})for(const auto& t:*messages)for(int i=0;i<2;++i)for(int j=0;j<2;++j)require(before[at++]==tcapi::get_elem(ctx,t,{i,j}),"difference must not mutate either input");
  for(std::size_t j=0;j<old.size();++j){next[j]=tcapi::copy(ctx,old[j]);tcapi::scale(ctx,next[j],E(7)*phase);}
  for(std::size_t j=0;j<ids.size();++j)for(int direction=0;direction<2;++direction){
    const auto edge=edges[ids[j]];
    if(owners.at(direction?edge.second:edge.first)==rank)
      tnbp::NormalizeMessage(ctx,next[j+direction*ids.size()],comm);
  }
  tnbp::BeliefPropagationDifference(ctx,edges,owners,old,next,ids,comm,residual);
  require(residual<tol,"scale / trace-phase invariance");
  std::vector<T> empty;std::vector<int> no_ids;std::vector<std::pair<int,int>> no_edges;
  tnbp::BeliefPropagationDifference(ctx,no_edges,owners,empty,empty,no_ids,comm,residual);
  require(residual==R(0),"empty graph residual");
  auto bad=std::vector<T>(old.size()+1);bool caught=false;
  try{tnbp::BeliefPropagationDifference(ctx,edges,owners,old,bad,ids,comm,residual);}catch(const std::invalid_argument&){caught=true;}
  require(caught,"message count validation");
  if(rank==0)std::cout<<"PASS successive difference type_bytes="<<sizeof(E)<<" ranks="<<size<<std::endl;
}
int main(int argc,char**argv){MPI_Init(&argc,&argv);try{
#ifdef TNBP_TEST_CUDA
 int rank,n=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);cudaGetDeviceCount(&n);require(n>0,"GPU required");cudaSetDevice(rank%n);
#endif
 run<double>(MPI_COMM_WORLD);run<std::complex<double>>(MPI_COMM_WORLD);
#ifndef TNBP_TEST_CUDA
 run<float>(MPI_COMM_WORLD);run<std::complex<float>>(MPI_COMM_WORLD);
#endif
 MPI_Finalize();return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,1);return 1;}}
