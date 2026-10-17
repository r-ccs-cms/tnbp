#include "tnbp/tnbp.h"
#include <complex>
#include <limits>
#include <sstream>
#include <iostream>
#ifdef TNBP_TEST_CUDA
#include <cuda_runtime.h>
template<class E> using Tensor = tcapi::cuda::Tensor<E>;
#else
template<class E> using Tensor = gqten::tensor<E>;
#endif
void check(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
template<class T> auto scalar(tnbp::context_handle_t<T>& ctx,const T& t) {
  return tcapi::get_elem(ctx,t,{});
}
template<class E> E phase(double x) {
  if constexpr(std::is_arithmetic_v<E>) return E(-1);
  else return std::polar(typename E::value_type(1),typename E::value_type(x));
}
template<class E> void run(MPI_Comm comm) {
  using T=Tensor<E>; using R=tnbp::real_t<T>;
  tnbp::context_handle_t<T> ctx; tcapi::create_context(ctx);
  const R tol=R(2000)*std::numeric_limits<R>::epsilon();
  std::stringstream fingerprint_input("hello");
  check(tnbp::CheckpointFingerprint(fingerprint_input)=="fnv1a64:a430d84680aabd0b","checkpoint fingerprint known vector");
  check(fingerprint_input.tellg()==0,"fingerprint preserves stream position");
  auto original=tcapi::eye<T>(ctx,2);
  tcapi::for_each_with_coors(ctx,original,[](auto& v,const auto& c){if(c[0]==c[1])v=E(c[0]+1);});
  for(double angle:{0.4,1.9,-2.4}) {
    auto m=tcapi::copy(ctx,original); tcapi::scale(ctx,m,E(3)*phase<E>(angle));
    T trace; tcapi::trace(ctx,m,{{0,1}},trace);
    check(std::abs(scalar(ctx,trace)-E(9)*phase<E>(angle))<tol*20,"logical trace/lazy scale");
    const auto norm=tcapi::norm(ctx,m);
    const E obs_before=tcapi::get_elem(ctx,m,{0,0})/scalar(ctx,trace);
    tnbp::FixMessageTracePhase(ctx,m);
    tcapi::trace(ctx,m,{{0,1}},trace);
    check(std::abs(scalar(ctx,trace)-E(9))<tol*20,"positive real trace");
    check(std::abs(tcapi::norm(ctx,m)-norm)<tol*norm,"norm preserved");
    check(std::abs(tcapi::get_elem(ctx,m,{0,0})/scalar(ctx,trace)-obs_before)<tol,"normalized observable invariant");
  }
  for(int invalid=0;invalid<4;++invalid) {
    auto m=tcapi::zeros<T>(ctx,{2,2});
    tcapi::for_each_with_coors(ctx,m,[&](auto& v,const auto& c){
      if(c[0]==c[1]) v=c[0]==0 ? E(1) : E(-1);
      if(invalid==1 && c[0]==1 && c[1]==1) v=E(-1+std::numeric_limits<R>::epsilon());
      if(invalid==2 && c[0]==0 && c[1]==0) v=E(std::numeric_limits<R>::infinity());
      if(invalid==3 && c[0]==0 && c[1]==0) v=E(std::numeric_limits<R>::quiet_NaN());
    });
    bool caught=false;try{tnbp::FixMessageTracePhase(ctx,m,"test-edge");}catch(const std::domain_error& e){caught=std::string(e.what()).find("test-edge")!=std::string::npos;}
    check(caught,"invalid trace diagnostic");
  }
  int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
  std::vector<std::pair<int,int>> edges={{0,1},{0,2},{1,2}};
  std::map<int,int> owners;std::vector<int> sites,edgeids;std::vector<T> v,e,f;
  for(int i=0;i<3;++i) {
    owners[i]=i%size;
    if(owners[i]==rank){sites.push_back(i); auto t=tcapi::fill<T>(ctx,{2,2,2},E(0));
      tcapi::for_each_with_coors(ctx,t,[i](auto& x,const auto& c){x=E(1+0.11*c[0]+0.23*c[1]+0.37*c[2]+0.07*i);});v.push_back(std::move(t));}
  }
  for(int k=0;k<3;++k)if(owners[edges[k].first]==rank || owners[edges[k].second]==rank)edgeids.push_back(k);
  for(int i=0;i<2*edgeids.size();++i)e.push_back(tcapi::eye<T>(ctx,2));
  R residual=0;
  for(int k=0;k<60;++k) {
    tnbp::BeliefPropagation(ctx,edges,v,sites,owners,edges,e,edgeids,comm,f);e=std::move(f);
    tnbp::BeliefPropagationCondition(ctx,edges,v,sites,owners,e,edgeids,comm,residual);
    if(residual<tol)break;
  }
  check(residual<tol,"small graph convergence");
  std::vector<T> observables;
  for(int i=0;i<3;++i) {
    auto op=tcapi::eye<T>(ctx,2);
    tcapi::for_each_with_coors(ctx,op,[](auto& x,const auto& c){if(c[0]==1 && c[1]==1)x=E(-1);});
    observables.push_back(std::move(op));
  }
  auto measured=tnbp::Measure(ctx,edges,v,sites,owners,e,edgeids,comm,std::vector<int>{0,1,2},observables);
  for(auto& m:e)tcapi::scale(ctx,m,phase<E>(1.9));
  auto rotated=tnbp::Measure(ctx,edges,v,sites,owners,e,edgeids,comm,std::vector<int>{0,1,2},observables);
  for(std::size_t i=0;i<measured.size();++i)check(std::abs(measured[i]-rotated[i])<tol,"BP observable phase invariance");
  R shifted=0;tnbp::BeliefPropagationCondition(ctx,edges,v,sites,owners,e,edgeids,comm,shifted);
  check(std::abs(shifted-residual)<tol,"phase invariant BP residual");
  std::stringstream checkpoint(std::ios::in|std::ios::out|std::ios::binary);
  tnbp::SaveTPS(ctx,checkpoint,v,sites,owners,e,edgeids);
  auto payload=checkpoint.str();
  for(bool legacy:{false,true}) {
    auto bytes=legacy ? payload.substr(0,payload.find("TNBP_BP_METADATA_V1")) : payload;
    std::stringstream input(bytes,std::ios::in|std::ios::binary);
    std::vector<T> vr,er,fr;std::vector<int> sr,ir;std::map<int,int> mr;
    tnbp::BPCheckpointMetadata metadata;
    tnbp::LoadTPS(ctx,input,vr,sr,mr,er,ir,&metadata);
    check(metadata.source_fingerprint.find("fnv1a64:")==0,"checkpoint source fingerprint");
    check(metadata.source_algorithm==(legacy?"legacy-unspecified":tnbp::BPAlgorithm),"algorithm provenance");
    auto resumed=tnbp::Measure(ctx,edges,vr,sr,mr,er,ir,comm,std::vector<int>{0,1,2},observables);
    for(std::size_t i=0;i<measured.size();++i)check(std::abs(measured[i]-resumed[i])<tol,"restart observable invariant");
    tnbp::BeliefPropagation(ctx,edges,v,sites,owners,edges,e,edgeids,comm,f);
    tnbp::BeliefPropagation(ctx,edges,vr,sr,mr,edges,er,ir,comm,fr);
    for(std::size_t j=0;j<f.size();++j)check(tcapi::close(ctx,f[j],fr[j],tol),"continuous/restart message match");
  }
  // Instantiate the public TPS overload too.
  tnbp::TensorProductState<T> state(ctx,v,edges,sites,owners,comm);
  tnbp::BeliefPropagation(ctx,state,edges,e,edgeids,f);
  if(rank==0)std::cout<<"PASS type_bytes="<<sizeof(E)<<" ranks="<<size<<" residual="<<residual<<std::endl;
}
int main(int argc,char**argv){MPI_Init(&argc,&argv);try{
#ifdef TNBP_TEST_CUDA
int rank;MPI_Comm_rank(MPI_COMM_WORLD,&rank);int n=0;cudaGetDeviceCount(&n);check(n>0,"CUDA device required");cudaSetDevice(rank%n);
#endif
run<double>(MPI_COMM_WORLD);run<std::complex<double>>(MPI_COMM_WORLD);
#ifndef TNBP_TEST_CUDA
run<float>(MPI_COMM_WORLD);run<std::complex<float>>(MPI_COMM_WORLD);
#endif
MPI_Finalize();return 0;}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,1);return 1;}}
