#include "tnbp/tnbp.h"
#include "tnbp/framework/timestamp.h"
#include <algorithm>
#include <complex>
#include <functional>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>
using C=std::complex<double>;
#ifdef TNBP_TEST_CUDA
using Tensor=tcapi::cuda::Tensor<C>;
#else
using Tensor=gqten::tensor<C>;
#endif
using Context=tnbp::context_handle_t<Tensor>;
using BondDim=tnbp::bond_dim_t<Tensor>;
using Real=tnbp::real_t<Tensor>;
using Elem=tnbp::elem_t<Tensor>;
const std::vector<std::pair<int,int>> edges{{0,1},{1,2},{3,4},{4,5}};
int checks=0;
void check(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);++checks;}
std::vector<Tensor> clone(Context& ctx,const std::vector<Tensor>& a){std::vector<Tensor> b;for(const auto&t:a)b.push_back(tcapi::copy(ctx,t));return b;}
std::string bytes(Context& ctx,const Tensor& t){std::ostringstream out(std::ios::binary);tcapi::save(ctx,t,out);return out.str();}
struct State {
  std::vector<Tensor> V,E;
  std::vector<int> sites,edge_idx;
  std::map<int,int> owners;
};
State clone(Context&ctx,const State&s){return {clone(ctx,s.V),clone(ctx,s.E),s.sites,s.edge_idx,s.owners};}
State initial(Context& ctx,MPI_Comm comm){
  int size;MPI_Comm_size(comm,&size);State s;
  if(size>=4)s.owners={{0,0},{1,1},{2,1},{3,2},{4,3},{5,3}};
  else if(size>=2)s.owners={{0,0},{1,1},{2,0},{3,0},{4,1},{5,0}};
  else for(int i=0;i<6;++i)s.owners.emplace(i,0);
  tnbp::InitTensorProductState(ctx,{0,1,2,3,4,5},edges,{2,2,2,2,2,2},s.V,s.sites,s.owners,s.E,s.edge_idx,comm);
  for(auto& v:s.V){auto shape=tcapi::shape(ctx,v);std::vector<Elem> d(2,1/std::sqrt(2.));v=tcapi::assign_from_range<Tensor>(ctx,shape,d.begin(),[](const auto&c){return c.back();});}
  return s;
}
void bp(Context&ctx,State&s,MPI_Comm comm){
  int rank;MPI_Comm_rank(comm,&rank);auto order=tnbp::MakeBPCommunicationOrder(edges,s.owners,rank);std::vector<Tensor> F;
  for(int i=0;i<50;++i){
    tnbp::UpdateBPMessages(ctx,edges,s.V,s.sites,s.E,s.edge_idx,comm,F);
    tnbp::ExchangeBPMessages(ctx,edges,s.owners,s.edge_idx,order,comm,F);
    Real residual=0;tnbp::BeliefPropagationDifference(ctx,edges,s.owners,s.E,F,s.edge_idx,comm,residual);s.E.swap(F);
    if(residual<1e-12)return;
  }
  throw std::runtime_error("batch fixture BP did not converge");
}
std::vector<C> normalized(std::vector<C> a){double n=0;for(auto z:a)n+=std::norm(z);n=std::sqrt(n);check(n>0 && std::isfinite(n),"valid state norm");for(auto&z:a)z/=n;return a;}
// Test-only collection/virtual enumeration, independent of BP measurement.
std::vector<C> dense(Context&ctx,const State&s,MPI_Comm comm){
  int rank;MPI_Comm_rank(comm,&rank);std::vector<std::vector<Elem>> data(6);std::vector<tnbp::shape_t<Tensor>> shape(6);std::vector<BondDim> dims(4,0);
  for(int site=0;site<6;++site){
    Tensor t;if(rank==s.owners.at(site)){auto i=std::find(s.sites.begin(),s.sites.end(),site)-s.sites.begin();t=tcapi::copy(ctx,s.V.at(i));}
    tnbp::MpiBcast(ctx,t,s.owners.at(site),comm);shape[site]=tcapi::shape(ctx,t);const auto& sh=shape[site];auto incident=tnbp::GetSurroundingBondIndex(site,edges);
    for(std::size_t k=0;k<incident.size();++k){auto& d=dims[incident[k]];check(!d || d==sh[k],"consistent virtual shape");d=sh[k];}
    data[site].resize(tcapi::size(ctx,t));tcapi::to_range(ctx,t,data[site].begin(),[sh](const auto&c){std::size_t a=0,k=1;for(std::size_t i=0;i<sh.size();++i){a+=c[i]*k;k*=sh[i];}return a;});
  }
  std::size_t assignments=1;for(auto d:dims)assignments*=d;std::vector<C> psi(64);
  for(std::size_t a=0;a<assignments;++a){auto rem=a;tnbp::elem_coors_t<Tensor> c(4);for(int e=0;e<4;++e){c[e]=rem%dims[e];rem/=dims[e];}
    for(int bits=0;bits<64;++bits){C value=1;for(int site=0;site<6;++site){std::size_t offset=0,stride=1;auto incident=tnbp::GetSurroundingBondIndex(site,edges);for(std::size_t k=0;k<incident.size();++k){offset+=c[incident[k]]*stride;stride*=shape[site][k];}offset+=((bits>>site)&1)*stride;value*=data[site][offset];}psi[bits]+=value;}
  }
  return normalized(psi);
}
double distance(const std::vector<C>&a,const std::vector<C>&b){C overlap=0;for(std::size_t i=0;i<a.size();++i)overlap+=std::conj(a[i])*b[i];const C phase=std::abs(overlap)>0 ? overlap/std::abs(overlap):C(1);double err=0;for(std::size_t i=0;i<a.size();++i)err+=std::norm(b[i]-phase*a[i]);return std::sqrt(err);}
std::vector<C> exact_step(const std::vector<C>&psi,const std::string&label,double x){std::vector<C> p(64);for(std::size_t bits=0;bits<64;++bits){auto out=bits;C v=psi[bits];for(int q=0;q<6;++q){auto op=label[5-q];bool bit=(bits>>q)&1;if(op=='X'||op=='Y')out^=std::size_t(1)<<q;if(op=='Y')v*=bit?C(0,-1):C(0,1);if(op=='Z'&&bit)v=-v;}p[out]+=v;}std::vector<C> r(64);for(int i=0;i<64;++i)r[i]=std::cosh(x)*psi[i]-std::sinh(x)*p[i];return normalized(r);}
auto batch(Context&ctx,State&s,const std::vector<int>&ids,MPI_Comm comm,std::vector<Real>*norms=nullptr){return tnbp::TruncationBatch(ctx,edges,ids,s.V,s.sites,s.owners,s.E,s.edge_idx,comm,100,1e-12,1e-12,norms);}
void all(Context&ctx,State&s,MPI_Comm comm,std::vector<BondDim>&dims,std::vector<Real>&err){tnbp::Truncation(ctx,edges,s.V,s.sites,s.owners,s.E,s.edge_idx,comm,100,1e-12,1e-12,dims,err);}
void equal_tensors(Context&ctx,const std::vector<Tensor>&a,const std::vector<Tensor>&b,const std::string&tag){check(a.size()==b.size(),tag+" count");for(std::size_t i=0;i<a.size();++i){check(tcapi::shape(ctx,a[i])==tcapi::shape(ctx,b[i]),tag+" shape");auto d=tcapi::linear_combine<Tensor>(ctx,{std::cref(a[i]),std::cref(b[i])},{C(1),C(-1)});check(tcapi::norm(ctx,d)<1e-11,tag+" values");}}
void api_contract(Context&ctx,MPI_Comm comm){
  int rank;MPI_Comm_rank(comm,&rank);auto s=initial(ctx,comm);
  auto p=tnbp::MakePauliEvolutionPlan(tnbp::NormalizePauliHamiltonian({{{.2,0},"IIIZYZ"},{{-.13,0},"YZYIII"}}),edges);
  auto t=tnbp::MakePauliLayerExponentialTPO<Tensor>(ctx,p,0,.04);
  tnbp::AbsorbTPO(ctx,t.tensors,p.graph.sites,edges,s.V,s.sites,s.owners,s.E,s.edge_idx,comm);bp(ctx,s,comm);
  auto before=clone(ctx,s);std::vector<Real> norms{123};auto empty=batch(ctx,s,{},comm,&norms);
  check(empty.edge_ids.empty()&&empty.bond_dims.empty()&&empty.truncation_errors.empty()&&norms==std::vector<Real>{123},"empty batch and norm no-op");
  equal_tensors(ctx,s.V,before.V,"empty V");equal_tensors(ctx,s.E,before.E,"empty E");
  for(const auto ids:std::vector<std::vector<int>>{{-1},{4},{0,0},{0,1}}){
    bool rejected=false;try{batch(ctx,s,ids,comm);}catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"invalid batch rejected before communication");
  }
  for(int mode=0;mode<3;++mode){
    bool rejected=false;
    try{tnbp::TruncationBatch(ctx,edges,{0},s.V,s.sites,s.owners,s.E,s.edge_idx,comm,
      mode==0 ? 0:100,mode==1 ? -1.0:1e-12,mode==2 ? std::numeric_limits<double>::quiet_NaN():1e-12);}
    catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"invalid numeric parameters rejected before communication");
  }
  auto missing_owner=s.owners;missing_owner.erase(0);bool rejected=false;
  try{tnbp::TruncationBatch(ctx,edges,{0},s.V,s.sites,missing_owner,s.E,s.edge_idx,comm,100,1e-12,1e-12);}
  catch(const std::invalid_argument&){rejected=true;}
  check(rejected,"missing selected endpoint owner rejected");
  bool involved=false;
  for(int e:{0,2}) involved=involved || s.owners.at(edges[e].first)==rank || s.owners.at(edges[e].second)==rank;
  for(int mode=0;mode<2;++mode){
    auto malformed=clone(ctx,s);
    if(mode==0 && !malformed.V.empty()) malformed.V.pop_back();
    if(mode==1 && !malformed.E.empty()) malformed.E.pop_back();
    bool rejected=false;
    try{batch(ctx,malformed,{0,2},comm);}catch(const std::invalid_argument&){rejected=true;}
    check(rejected==involved,"malformed tensor counts rejected locally on participants");
  }
  for(std::size_t i=0;i<s.V.size();++i)check(bytes(ctx,s.V[i])==bytes(ctx,before.V[i]),"invalid calls preserve V bytes");
  for(std::size_t i=0;i<s.E.size();++i)check(bytes(ctx,s.E[i])==bytes(ctx,before.E[i]),"invalid calls preserve E bytes");
  auto result=batch(ctx,s,{2,0},comm,&norms);std::vector<int> expected;
  for(int e:{0,2})if(s.owners.at(edges[e].first)==rank||s.owners.at(edges[e].second)==rank)expected.push_back(e);
  check(result.edge_ids==expected,"local selected result IDs canonical order");check(result.bond_dims.size()==expected.size()&&result.truncation_errors.size()==expected.size(),"selected-only result sizes");
  std::set<int> touched{0,1,3,4};
  for(std::size_t i=0;i<s.sites.size();++i)if(!touched.count(s.sites[i]))check(bytes(ctx,s.V[i])==bytes(ctx,before.V[i]),"unselected site bytes unchanged");
  for(std::size_t i=0;i<s.edge_idx.size();++i)if(s.edge_idx[i]==1||s.edge_idx[i]==3){check(bytes(ctx,s.E[i])==bytes(ctx,before.E[i]),"unselected forward E unchanged");check(bytes(ctx,s.E[i+s.edge_idx.size()])==bytes(ctx,before.E[i+s.edge_idx.size()]),"unselected reverse E unchanged");}
  check(expected.empty() ? norms==std::vector<Real>{123}:norms.size()>1,"norms appended only on participants");
  // Owners of the other component deliberately omit this API call. If a
  // collective is accidentally introduced, the following test Bcasts hang.
  auto local_only=clone(ctx,before);if(rank==s.owners.at(0)||rank==s.owners.at(1))batch(ctx,local_only,{0},comm);
  auto psi=dense(ctx,local_only,comm);check(psi.size()==64,"nonparticipants need not enter batch");
  auto whole=clone(ctx,before),split=clone(ctx,before);std::vector<BondDim>dims;std::vector<Real>errors;all(ctx,whole,comm,dims,errors);
  std::vector<BondDim> split_dims(s.edge_idx.size());std::vector<Real>split_errors(s.edge_idx.size());
  for(const auto ids:{std::vector<int>{2,0},std::vector<int>{3,1}}){auto r=batch(ctx,split,ids,comm);for(std::size_t i=0;i<r.edge_ids.size();++i){auto pos=std::find(s.edge_idx.begin(),s.edge_idx.end(),r.edge_ids[i])-s.edge_idx.begin();split_dims[pos]=r.bond_dims[i];split_errors[pos]=r.truncation_errors[i];}}
  check(dims==split_dims,"batch sweep retained dimensions equal all-edge call");for(std::size_t i=0;i<errors.size();++i)check(std::abs(errors[i]-split_errors[i])<1e-14,"batch sweep errors equal all-edge call");
  equal_tensors(ctx,whole.V,split.V,"full versus batch V");equal_tensors(ctx,whole.E,split.E,"full versus batch E");
}
void real_tensors(Context& complex_ctx,MPI_Comm comm){
  using RealTensor=tcapi::real_ten_t<Tensor>;
  tcapi::context_handle_t<RealTensor> ctx;tcapi::create_context(ctx);
  auto template_state=initial(complex_ctx,comm);
  std::vector<RealTensor> V,E;std::vector<int> sites,edge_idx;
  tnbp::InitTensorProductState(ctx,{0,1,2,3,4,5},edges,{2,2,2,2,2,2},V,sites,template_state.owners,E,edge_idx,comm);
  auto r=tnbp::TruncationBatch(ctx,edges,{0,2},V,sites,template_state.owners,E,edge_idx,comm,100,1e-12,1e-12);
  for(auto d:r.bond_dims)check(d==1,"real tensor retained dimension");
  for(auto e:r.truncation_errors)check(e==0,"real tensor truncation error");
}
double evolution(Context&ctx,MPI_Comm comm,const std::string&first,const std::string&second,bool refresh){
  auto s=initial(ctx,comm);auto p=tnbp::MakePauliEvolutionPlan(tnbp::NormalizePauliHamiltonian({{{.2,0},first},{{-.13,0},second}}),edges);
  check(p.layers.size()==1,"independent three-site operators in one layer");auto t=tnbp::MakePauliLayerExponentialTPO<Tensor>(ctx,p,0,.04,s.sites);
  check(t.tensors.size()==s.sites.size() && std::abs(t.log_scale-.0132)<1e-15,"local TPO count and whole-layer scale, including empty ranks");
  auto reference=dense(ctx,s,comm);double error=0;
  for(int step=0;step<10;++step){
    reference=exact_step(reference,first,.008);reference=exact_step(reference,second,-.0052);
    tnbp::AttachTPO(ctx,edges,t.tensors,s.V,s.sites,s.owners,s.E,s.edge_idx,comm);bp(ctx,s,comm);
    auto batches=tnbp::MakePauliLayerTruncationBatches(p,0);
    check(batches==std::vector<std::vector<int>>{{1,3},{0,2}},"leaf-side automatic batches");
    for(std::size_t b=0;b<batches.size();++b){if(b && refresh)bp(ctx,s,comm);batch(ctx,s,batches[b],comm);}
    const double current=distance(reference,dense(ctx,s,comm));error=std::max(error,current);
    if(refresh)check(current<2e-10,"multiple three-site operations: batch/BP/batch versus exact state");
  }
  return error;
}
int main(int argc,char**argv){
  MPI_Init(&argc,&argv);int rank;MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  try{
#ifdef TNBP_TEST_CUDA
    int count=0;cudaGetDeviceCount(&count);check(count==1,"one visible GPU per rank");cudaSetDevice(0);
#endif
    Context ctx;tcapi::create_context(ctx);api_contract(ctx,MPI_COMM_WORLD);real_tensors(ctx,MPI_COMM_WORLD);
    double exact=evolution(ctx,MPI_COMM_WORLD,"IIIZYZ","YZYIII",true);
    double long_range=evolution(ctx,MPI_COMM_WORLD,"IIIZIZ","ZIYIII",true);
    double frozen=evolution(ctx,MPI_COMM_WORLD,"IIIZIZ","ZIYIII",false);
    if(rank==0)std::cout<<tnbp::make_timestamp()<<" # truncation_batch PASS checks="<<checks<<" three_site_error="<<exact<<" long_range_error="<<long_range<<" frozen_environment_error="<<frozen<<std::endl;
  }catch(const std::exception&e){std::cerr<<tnbp::make_timestamp()<<" rank="<<rank<<" "<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,1);}
  MPI_Finalize();
}
