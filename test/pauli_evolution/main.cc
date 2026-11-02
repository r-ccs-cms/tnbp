#include "tnbp/operator/pauli_exponential.h"
#include "tnbp/ptns/function/apply.h"
#include "tnbp/ptns/function/init.h"
#include "tnbp/framework/timestamp.h"
#include <iostream>
#include <functional>
#include <mpi.h>
using C=std::complex<double>;
using Dense=std::vector<C>;
using Edges=std::vector<std::pair<int,int>>;
#ifdef TNBP_TEST_CUDA
using Tensor=tcapi::cuda::Tensor<C>;
#else
using Tensor=gqten::tensor<C>;
#endif
using Context=tnbp::context_handle_t<Tensor>;
int checks=0;
double max_operator_error=0;
std::array<double,2> convergence_ratios{};
void check(bool ok,const std::string& message) {
  if(!ok) throw std::runtime_error(message);
  ++checks;
}
void rejects(const std::function<void()>& f,const std::string& message) {
  bool rejected=false;
  try {f();} catch(const std::invalid_argument&) {rejected=true;} catch(const std::overflow_error&) {rejected=true;}
  check(rejected,message);
}
Dense identity(std::size_t n) {Dense a(n*n);for(std::size_t i=0;i<n;++i)a[i*n+i]=1;return a;}
Dense mul(const Dense& a,const Dense& b,std::size_t n) {
  Dense c(n*n);
  for(std::size_t i=0;i<n;++i)for(std::size_t k=0;k<n;++k)for(std::size_t j=0;j<n;++j)c[i*n+j]+=a[i*n+k]*b[k*n+j];
  return c;
}
double distance(const Dense& a,const Dense& b) {
  if(a.size()!=b.size()) throw std::runtime_error("dense shape mismatch");
  double d=0;for(std::size_t i=0;i<a.size();++i)d=std::max(d,std::abs(a[i]-b[i]));return d;
}
// Independent bit-basis action: no TPO local matrix helper or tensor contractions.
Dense pauli_matrix(const std::string& label) {
  const std::size_t n=std::size_t(1)<<label.size(); Dense p(n*n);
  for(std::size_t input=0;input<n;++input) {
    auto output=input; C amplitude=1;
    for(std::size_t q=0;q<label.size();++q) {
      const bool bit=(input>>q)&1;
      switch(label[label.size()-1-q]) {
        case 'I':break;
        case 'X':output^=std::size_t(1)<<q;break;
        case 'Y':output^=std::size_t(1)<<q;amplitude*=bit ? C(0,-1) : C(0,1);break;
        case 'Z':if(bit)amplitude=-amplitude;break;
        default:throw std::runtime_error("reference label");
      }
    }
    p[output*n+input]=amplitude;
  }
  return p;
}
Dense exponential(const std::string& label,double x,bool scaled=false) {
  auto p=pauli_matrix(label);const auto n=std::size_t(1)<<label.size();auto a=identity(n);
  double c=std::cosh(x),s=-std::sinh(x);
  if(scaled) {const double u=std::exp(-2*std::abs(x));c=(1+u)/2;s=-std::copysign((1-u)/2,x);}
  for(std::size_t i=0;i<a.size();++i)a[i]=c*a[i]+s*p[i];return a;
}
Dense tpo_dense(Context& ctx,const tnbp::PauliEvolutionPlan& p,const tnbp::ScaledPauliTPO<Tensor>& tpo) {
  const auto n=std::size_t(1)<<p.hamiltonian.num_sites;
  std::vector<std::vector<C>> data;std::vector<tnbp::shape_t<Tensor>> shapes;
  std::vector<std::size_t> dims(p.graph.edges.size(),0);
  for(std::size_t v=0;v<tpo.tensors.size();++v) {
    const auto shape=tcapi::shape(ctx,tpo.tensors[v]);shapes.push_back(shape);
    check(shape.size()==p.graph.incident_edge_ids[v].size()+2,"TPO rank");
    check(shape[shape.size()-1]==2 && shape[shape.size()-2]==2,"TPO physical shape");
    for(std::size_t k=0;k+2<shape.size();++k) {
      const auto e=p.graph.incident_edge_ids[v][k];
      check(!dims[e] || dims[e]==static_cast<std::size_t>(shape[k]),"TPO matching edge dimensions");
      dims[e]=shape[k];check(dims[e]<=2,"TPO edge cap 2");
    }
    std::vector<C> flat(tcapi::size(ctx,tpo.tensors[v]));
    tcapi::to_range(ctx,tpo.tensors[v],flat.begin(),[shape](const auto& c){
      std::size_t offset=0,stride=1;for(std::size_t k=0;k<shape.size();++k){offset+=c[k]*stride;stride*=shape[k];}return offset;
    }); data.push_back(std::move(flat));
  }
  std::size_t assignments=1;for(auto d:dims)assignments*=d;
  Dense a(n*n);
  for(std::size_t b=0;b<assignments;++b) {
    auto rem=b;std::vector<std::size_t> coord(dims.size());
    for(std::size_t e=0;e<dims.size();++e){coord[e]=rem%dims[e];rem/=dims[e];}
    for(std::size_t o=0;o<n;++o)for(std::size_t i=0;i<n;++i) {
      C product=1;
      for(std::size_t v=0;v<data.size();++v) {
        std::size_t offset=0,stride=1;
        for(std::size_t k=0;k<p.graph.incident_edge_ids[v].size();++k) {
          offset+=coord[p.graph.incident_edge_ids[v][k]]*stride;stride*=shapes[v][k];
        }
        offset+=((i>>v)&1)*stride;stride*=2;offset+=((o>>v)&1)*stride;
        product*=data[v][offset];
      }
      a[o*n+i]+=product;
    }
  }
  return a;
}
tnbp::PauliEvolutionPlan plan(std::vector<pauli::Term<C>> terms,const Edges& edges) {
  return tnbp::MakePauliEvolutionPlan(tnbp::NormalizePauliHamiltonian(terms),edges);
}
// Independent contract checks: removing an edge identifies all edges below it.
// This does not reproduce the planner's parent/ready-queue algorithm.
void batch_contract(const tnbp::PauliEvolutionPlan& p) {
  for(std::size_t l=0;l<p.layers.size();++l) {
    std::set<int> layer_sites;
    for(auto t:p.layers[l].term_ids)for(int v:p.routes[t].tree_sites)
      check(layer_sites.insert(v).second,"layer trees share no sites, including routing sites");
    auto batches=tnbp::MakePauliLayerTruncationBatches(p,l);
    std::map<std::size_t,std::size_t> position;
    for(std::size_t b=0;b<batches.size();++b) {
      std::set<int> occupied;
      for(int e:batches[b]) {
        auto [a,c]=p.graph.edges.at(e);
        check(occupied.insert(a).second && occupied.insert(c).second,"batch has no site conflict");
        check(position.emplace(e,b).second,"edge compressed exactly once");
      }
      for(auto t:p.layers[l].term_ids) {
        std::size_t count=0;
        for(auto e:p.routes[t].tree_edge_ids)count+=std::find(batches[b].begin(),batches[b].end(),static_cast<int>(e))!=batches[b].end();
        check(count<=1,"at most one edge of each tree per batch");
      }
    }
    check(position.size()==p.layers[l].edge_ids.size(),"all layer edges covered");
    for(auto t:p.layers[l].term_ids) {
      const auto& route=p.routes[t];
      for(auto cut:route.tree_edge_ids) {
        std::set<int> root_side{route.root_site};
        bool changed=true;
        while(changed) {
          changed=false;
          for(auto e:route.tree_edge_ids)if(e!=cut) {
            auto [a,b]=p.graph.edges[e];
            if(root_side.count(a))changed=root_side.insert(b).second || changed;
            if(root_side.count(b))changed=root_side.insert(a).second || changed;
          }
        }
        for(auto below:route.tree_edge_ids)if(below!=cut) {
          auto [a,b]=p.graph.edges[below];
          if(!root_side.count(a) && !root_side.count(b))
            check(position.at(below)<position.at(cut),"descendant edge finishes in earlier batch than ancestor");
        }
      }
    }
  }
}
void planning() {
  auto h=tnbp::NormalizePauliHamiltonian({{{1,2}," xi "},{{2,-2},"XI"},{{3,0},"II"},{{0,0},"ZI"}});
  check(h.num_sites==2 && h.terms.size()==1 && h.terms[0].coefficient==3 && h.identity_coefficient==3,"normalization merge/identity");
  rejects([]{tnbp::NormalizePauliHamiltonian({{{1,1},"Y"}});},"reject nonreal");
  rejects([]{tnbp::NormalizePauliHamiltonian({{{1,0},"A"}});},"reject label");
  rejects([]{tnbp::NormalizePauliHamiltonian({{{1,0},"X"},{{1,0},"XX"}});},"reject length");
  rejects([]{tnbp::NormalizePauliHamiltonian({{{INFINITY,0},"X"}});},"reject infinity");
  rejects([]{tnbp::NormalizePauliHamiltonian({});},"reject empty");
  rejects([]{plan({{{1,0},"XX"}},{{0,1},{1,0}});},"reject reverse duplicate");
  rejects([]{plan({{{1,0},"XX"}},{{0,0}});},"reject self edge");
  rejects([]{plan({{{1,0},"XX"}},{{0,2}});},"reject out of range");
  rejects([]{plan({{{1,0},"XIX"}},{});},"reject disconnected support");
  auto empty=plan({{{0,0},"IIII"}},{});
  check(empty.graph.sites.size()==4 && empty.layers.empty() && tnbp::MakeTrotterActions(empty).empty(),"zero Hamiltonian retains sites");
  auto unit=plan({{{2,0},"I"}},{});check(unit.layers.empty() && unit.hamiltonian.identity_coefficient==2,"identity has no layer");
  auto isolated=plan({{{1,0},"YII"}},{{0,1}});check(isolated.routes[0].tree_sites==std::vector<int>{2},"isolated local term");
  auto tied=plan({{{1,0},"ZIIZ"}},{{2,3},{0,2},{1,3},{0,1}});
  check(tied.routes[0].tree_sites==std::vector<int>({0,1,3}),"shortest path tie uses site order");
  check(tied.graph.edges==Edges({{2,3},{0,2},{1,3},{0,1}}),"preserve edge order without sorting");
  auto again=plan({{{1,0},"ZIIZ"}},{{2,3},{0,2},{1,3},{0,1}});
  check(tied.plan_digest==again.plan_digest && tied.plan_digest.size()==64,"stable digest");
  auto branched=plan({{{1,0},"ZZZI"}},{{0,1},{0,2},{0,3}});
  check(branched.routes[0].tree_edge_ids.size()==3 && branched.routes[0].tree_sites.size()==4,"branched support tree");
  auto shared=plan({{{1,0},"XXIX"},{{1,0},"ZZII"}},{{0,1},{1,2},{1,3},{2,3}});
  check(shared.layers.size()==2,"site-sharing commuting terms split despite disjoint edges");
  auto noncommute=plan({{{1,0},"X"},{{1,0},"Z"}},{});check(noncommute.layers.size()==2,"noncommuting single sites split");
  auto overlap=plan({{{1,0},"ZIZ"},{{1,0},"IZZ"}},{{0,1},{1,2}});check(overlap.layers.size()==2,"commuting overlapping routes split");
  auto leaf=plan({{{1,0},"ZZZZZ"}},{{0,1},{0,2},{1,3},{1,4}});
  check(tnbp::MakePauliLayerTruncationBatches(leaf,0)==std::vector<std::vector<int>>{{1},{2},{3},{0}},"branch: child edges before parent, smallest ready ID first");
  auto forest=plan({{{.17,0},"IIIZYIZ"},{{-.11,0},"YZYIIII"}},{{0,1},{1,2},{1,3},{4,5},{5,6}});
  check(tnbp::MakePauliLayerTruncationBatches(forest,0)==std::vector<std::vector<int>>{{1,4},{2,3},{0}},"different trees overlap one edge each");
  auto conflict=plan({{{1,0},"IIXXI"},{{1,0},"ZZIII"}},{{0,1},{0,2},{0,3},{0,4}});
  check(conflict.layers.size()==2,"intermediate identity site shared by disjoint supports splits layers");
  check(tnbp::MakePauliLayerTruncationBatches(conflict,0)==std::vector<std::vector<int>>{{1},{0}},"first disjoint-support tree remains leaf-first");
  check(tnbp::MakePauliLayerTruncationBatches(conflict,1)==std::vector<std::vector<int>>{{3},{2}},"second tree is in a separate layer");
  auto onsite_route=plan({{{1,0},"IXI"},{{1,0},"ZIZ"}},{{0,1},{1,2}});
  check(onsite_route.layers.size()==2,"onsite term cannot share an intermediate route site");
  auto common_z=plan({{{1,0},"IXZ"},{{1,0},"YIZ"}},{{0,1},{0,2}});
  check(common_z.layers.size()==2,"commuting common-site Pauli terms use separate layers");
  check(common_z.plan_version==2,"site-disjoint layer semantics use plan version 2");
  check(tnbp::MakePauliLayerTruncationBatches(isolated,0).empty(),"single-site term has no compression batches");
  rejects([&]{tnbp::MakePauliLayerTruncationBatches(leaf,1);},"reject invalid batch layer");
  batch_contract(forest);batch_contract(conflict);batch_contract(isolated);
  Edges permuted{{0,1},{0,2},{1,3},{1,4}};
  do{batch_contract(plan({{{1,0},"ZZZZZ"}},permuted));}while(std::next_permutation(permuted.begin(),permuted.end()));
  auto one=tnbp::MakeTrotterActions(isolated);check(one.size()==1 && one[0].time_fraction==1,"one layer full step");
  auto acts=tnbp::MakeTrotterActions(noncommute);check(acts.size()==3 && acts[0].layer_id==0 && acts[0].time_fraction==0.5 && acts[1].layer_id==1 && acts[1].time_fraction==1 && acts[2].layer_id==0,"second order sequence");
  rejects([&]{tnbp::MakeTrotterActions(noncommute,3);},"reject Trotter order");
}
void operators(Context& ctx) {
  const std::vector<std::pair<std::string,Edges>> cases={
    {"X",{}},{"Y",{}},{"Z",{}},{"YII",{{0,1}}},{"ZY",{{0,1}}},
    {"ZIIZ",{{2,3},{0,1},{1,2}}},{"YIXZ",{{0,1},{1,2},{2,3}}},
    {"XYZI",{{0,1},{0,2},{0,3}}},{"ZZZYZ",{{0,1},{1,2},{1,3},{1,4}}},
    {"YIYIY",{{0,1},{1,2},{2,3},{3,4},{4,0}}}};
  for(const auto& [label,edges]:cases) for(double c:{-0.7,0.3}) for(double dt:{0.0,1e-14,0.17,0.8}) {
    const auto p=plan({{{c,0},label}},edges); auto t=tnbp::MakePauliTermExponentialTPO<Tensor>(ctx,p,0,dt);
    auto actual=tpo_dense(ctx,p,t);auto expected=exponential(label,c*dt,true);
    const double error=distance(actual,expected);
    max_operator_error=std::max(max_operator_error,error);
    check(error<2e-13,"single exponential "+label);
    check(std::abs(t.log_scale-std::abs(c*dt))<1e-15,"single scale");
  }
  auto tiny=plan({{{1,0},"X"}},{});auto t=tnbp::MakePauliTermExponentialTPO<Tensor>(ctx,tiny,0,1e-20);
  auto a=tpo_dense(ctx,tiny,t);check(std::abs(a[1]+1e-20)<1e-34,"tiny coefficient avoids cancellation");
  auto tiny_route=plan({{{1,0},"XIZ"}},{{0,1},{1,2}});
  auto tiny_tpo=tnbp::MakePauliTermExponentialTPO<Tensor>(ctx,tiny_route,0,1e-300);
  auto tiny_matrix=tpo_dense(ctx,tiny_route,tiny_tpo);
  check(std::abs(tiny_matrix[4*8]+1e-300)<1e-313,"distributed tiny branch remains nonzero");
  for(double c:{-1000.,1000.}) {
    auto p=plan({{{c,0},"YIY"}},{{0,1},{1,2}});auto t=tnbp::MakePauliTermExponentialTPO<Tensor>(ctx,p,0,1);
    check(distance(tpo_dense(ctx,p,t),exponential("YIY",c,true))<1e-14,"large scaled exponential");
  }
  auto p=plan({{{0.4,0},"XXIX"},{{-0.7,0},"ZZII"},{{0.2,0},"IIXI"}},{{0,1},{1,2},{1,3},{2,3}});
  for(std::size_t l=0;l<p.layers.size();++l) {
    auto t=tnbp::MakePauliLayerExponentialTPO<Tensor>(ctx,p,l,0.2);
    auto expected=identity(16);double scale=0;
    for(auto id:p.layers[l].term_ids) {
      const auto& term=p.hamiltonian.terms[id];scale+=std::abs(0.2*term.coefficient);
      expected=mul(exponential(term.label,0.2*term.coefficient,true),expected,16);
    }
    const double error=distance(tpo_dense(ctx,p,t),expected);
    max_operator_error=std::max(max_operator_error,error);
    check(error<2e-13,"site-disjoint layer ordered product");
    check(std::abs(t.log_scale-scale)<1e-15,"layer scale sum");
  }
  rejects([&]{tnbp::MakePauliTermExponentialTPO<Tensor>(ctx,p,999,1);},"invalid term ID");
  rejects([&]{tnbp::MakePauliLayerExponentialTPO<Tensor>(ctx,p,999,1);},"invalid layer ID");
  rejects([&]{tnbp::MakePauliTermExponentialTPO<Tensor>(ctx,p,0,-1);},"negative dt");
  rejects([&]{tnbp::MakePauliTermExponentialTPO<Tensor>(ctx,p,0,INFINITY);},"nonfinite dt");
  auto huge=plan({{{1e308,0},"X"}},{});
  rejects([&]{tnbp::MakePauliTermExponentialTPO<Tensor>(ctx,huge,0,10);},"x overflow");
}
// Inspect the required factorization, which a contracted-matrix test alone
// cannot distinguish from the former root-only coefficient placement.
template<class TenT>
void support_weights(tnbp::context_handle_t<TenT>& ctx) {
  using Real=tnbp::real_t<TenT>;
  const double tolerance=64*std::numeric_limits<Real>::epsilon();
  const std::vector<std::pair<std::string,Edges>> cases={
    {"ZIIZ",{{0,1},{1,2},{2,3}}},
    {"XYZI",{{0,1},{0,2},{0,3}}},
    {"ZZZYZ",{{0,1},{1,2},{1,3},{1,4}}},
    {"XIZIY",{{0,1},{1,2},{2,3},{3,4}}}};
  for(const auto& [label,edges]:cases)for(double coefficient:{-.7,.3}) {
    auto p=plan({{{coefficient,0},label}},edges);
    const auto& route=p.routes[0];const double dt=.25;
    auto t=tnbp::MakePauliTermExponentialTPO<TenT>(ctx,p,0,dt);
    std::array<C,2> products{C(1),C(1)};
    std::array<double,2> reference_magnitudes{};bool first=true;
    for(const auto& local:route.local_sites) {
      const auto& tensor=t.tensors[local.site];auto shape=tcapi::shape(ctx,tensor);
      std::array<C,2> amplitudes{};
      for(tnbp::elem_coor_t<TenT> branch=0;branch<2;++branch) {
        auto reference=pauli_matrix(std::string(1,branch ? local.op:'I'));
        tnbp::elem_coors_t<TenT> coor(shape.size(),0);
        for(auto leg:local.tree_legs)coor[leg]=branch;
        for(tnbp::elem_coor_t<TenT> in=0;in<2;++in)for(tnbp::elem_coor_t<TenT> out=0;out<2;++out) {
          coor[shape.size()-2]=in;coor[shape.size()-1]=out;
          amplitudes[branch]+=std::conj(reference[out*2+in])*C(tcapi::get_elem(ctx,tensor,coor))*.5;
        }
        for(tnbp::elem_coor_t<TenT> in=0;in<2;++in)for(tnbp::elem_coor_t<TenT> out=0;out<2;++out) {
          coor[shape.size()-2]=in;coor[shape.size()-1]=out;
          check(std::abs(C(tcapi::get_elem(ctx,tensor,coor))-amplitudes[branch]*reference[out*2+in])<tolerance,
                "local branch is its Pauli matrix times one site weight");
        }
      }
      if(local.op=='I') {
        check(std::abs(amplitudes[0]-C(1))<tolerance && std::abs(amplitudes[1]-C(1))<tolerance,
              "identity routing and branch sites have no weight");
      } else {
        for(std::size_t branch=0;branch<2;++branch) {
          if(first)reference_magnitudes[branch]=std::abs(amplitudes[branch]);
          check(std::abs(std::abs(amplitudes[branch])-reference_magnitudes[branch])<tolerance,
                "all Pauli support sites have equal branch magnitudes");
          products[branch]*=amplitudes[branch];
        }
        check(amplitudes[0].real()>0 && std::abs(amplitudes[0].imag())<tolerance,"positive identity-branch weight");
        const bool negative=coefficient>0 && local.site==route.root_site;
        check((amplitudes[1].real()<0)==negative && std::abs(amplitudes[1].imag())<tolerance,
              "negative branch sign occurs only at the deterministic Pauli root");
        first=false;
      }
      if(local.tree_legs.size()>1) {
        tnbp::elem_coors_t<TenT> mismatch(shape.size(),0);mismatch[local.tree_legs.back()]=1;
        check(std::abs(tcapi::get_elem(ctx,tensor,mismatch))==0,"branch point synchronizes virtual flags");
      }
    }
    const double x=dt*coefficient,scale=std::exp(-std::abs(x));
    check(std::abs(products[0]-C(scale*std::cosh(x)))<tolerance,"product of site weights gives scaled cosh");
    check(std::abs(products[1]+C(scale*std::sinh(x)))<tolerance,"product of signed site weights gives scaled minus-sinh");
  }
}
// Taylor exp of a complete dense noncommuting Hamiltonian, independently of Trotter.
Dense exact_exp(Dense h,std::size_t n,double time) {
  auto sum=identity(n),term=identity(n);
  for(auto& x:h)x*=-time;
  for(int k=1;k<=50;++k) {term=mul(h,term,n);for(auto& x:term)x/=double(k);for(std::size_t i=0;i<sum.size();++i)sum[i]+=term[i];}
  return sum;
}
void trotter(Context& ctx) {
  auto p=plan({{{0.7,0},"XI"},{{-0.4,0},"ZZ"},{{0.2,0},"IY"},{{0.13,0},"II"}},{{0,1}});
  Dense h(16);for(const auto& term:p.hamiltonian.terms) {auto a=pauli_matrix(term.label);for(std::size_t i=0;i<h.size();++i)h[i]+=term.coefficient*a[i];}
  for(int i=0;i<4;++i)h[4*i+i]+=p.hamiltonian.identity_coefficient;
  auto exact=exact_exp(h,4,0.4);
  for(int order:{1,2}) {
    std::vector<double> errors;
    for(int steps:{4,8,16}) {
      double dt=0.4/steps;auto single=identity(4);double log_scale=-dt*p.hamiltonian.identity_coefficient;
      for(auto action:tnbp::MakeTrotterActions(p,order)) {
        auto t=tnbp::MakePauliLayerExponentialTPO<Tensor>(ctx,p,action.layer_id,dt*action.time_fraction);
        single=mul(tpo_dense(ctx,p,t),single,4);log_scale+=t.log_scale;
      }
      for(auto& x:single)x*=std::exp(log_scale);
      auto full=identity(4);for(int k=0;k<steps;++k)full=mul(single,full,4);
      errors.push_back(distance(full,exact));
    }
    const double ratio=errors[1]/errors[2];
    convergence_ratios[order-1]=ratio;
    check(ratio>(order==1 ? 1.8:3.7) && ratio<(order==1 ? 2.2:4.3),"Trotter convergence order "+std::to_string(order));
  }
}
std::string tensor_bytes(Context& ctx,const Tensor& tensor) {
  std::ostringstream out(std::ios::binary);tcapi::save(ctx,tensor,out);return out.str();
}
void local_generation(Context& ctx) {
  for(const auto& p:std::vector<tnbp::PauliEvolutionPlan>{
      plan({{{.2,0},"IZYIZ"}},{{0,1},{1,2},{1,3}}),
      plan({{{.2,0},"IIIZYIZ"},{{-.13,0},"YZYIIII"}},{{0,1},{1,2},{1,3},{4,5},{5,6}})}) {
    const std::vector<int> subset=p.graph.sites.size()==5 ? std::vector<int>{4,3,1}:std::vector<int>{6,2,5};
    for(double dt:{0.,1e-20,.17,1000.}) {
      auto check_subset=[&](const auto& all,const auto& local,const auto& empty) {
        check(all.tensors.size()==p.graph.sites.size() && local.tensors.size()==subset.size(),"only requested TPO tensors generated");
        check(empty.tensors.empty(),"empty local site list creates no tensors");
        check(local.log_scale==all.log_scale && empty.log_scale==all.log_scale,"whole-operator scalar independent of owned sites");
        for(std::size_t i=0;i<subset.size();++i)
          check(tensor_bytes(ctx,local.tensors[i])==tensor_bytes(ctx,all.tensors[subset[i]]),"local tensor equals full tensor byte-for-byte in requested order");
      };
      for(std::size_t t=0;t<p.hamiltonian.terms.size();++t)
        check_subset(tnbp::MakePauliTermExponentialTPO<Tensor>(ctx,p,t,dt),
          tnbp::MakePauliTermExponentialTPO<Tensor>(ctx,p,t,dt,subset),
          tnbp::MakePauliTermExponentialTPO<Tensor>(ctx,p,t,dt,{}));
      for(std::size_t l=0;l<p.layers.size();++l)
        check_subset(tnbp::MakePauliLayerExponentialTPO<Tensor>(ctx,p,l,dt),
          tnbp::MakePauliLayerExponentialTPO<Tensor>(ctx,p,l,dt,subset),
          tnbp::MakePauliLayerExponentialTPO<Tensor>(ctx,p,l,dt,{}));
    }
    for(const auto& sites:std::vector<std::vector<int>>{{-1},{static_cast<int>(p.graph.sites.size())},{1,1}}) {
      rejects([&]{tnbp::MakePauliTermExponentialTPO<Tensor>(ctx,p,0,.1,sites);},"invalid local term site list");
      rejects([&]{tnbp::MakePauliLayerExponentialTPO<Tensor>(ctx,p,0,.1,sites);},"invalid local layer site list");
    }
  }
}
void single_precision() {
#ifdef TNBP_TEST_CUDA
  using FTensor=tcapi::cuda::Tensor<std::complex<float>>;
#else
  using FTensor=gqten::tensor<std::complex<float>>;
#endif
  tnbp::context_handle_t<FTensor> ctx;tcapi::create_context(ctx);
  support_weights<FTensor>(ctx);
  auto p=plan({{{-0.7,0},"Y"}},{});
  auto t=tnbp::MakePauliLayerExponentialTPO<FTensor>(ctx,p,0,0.3);
  auto local=tnbp::MakePauliLayerExponentialTPO<FTensor>(ctx,p,0,.3,{0});
  auto empty=tnbp::MakePauliLayerExponentialTPO<FTensor>(ctx,p,0,.3,{});
  check(local.tensors.size()==1 && empty.tensors.empty() && local.log_scale==t.log_scale && empty.log_scale==t.log_scale,"complex<float> local generation and common scale");
  const auto expected=exponential("Y",-0.21,true);
  for(tnbp::elem_coor_t<FTensor> i=0;i<2;++i)for(tnbp::elem_coor_t<FTensor> o=0;o<2;++o) {
    const auto actual=tcapi::get_elem(ctx,t.tensors[0],{i,o});
    check(std::abs(C(actual)-expected[o*2+i])<1e-7,"complex<float> exponential");
    check(tcapi::get_elem(ctx,local.tensors[0],{i,o})==actual,"complex<float> local/global values identical");
  }
}
void absorption(Context& ctx,MPI_Comm comm) {
  int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
  auto p=plan({{{0.37,0},"YZ"}},{{0,1}});
  std::map<int,int> owners{{0,0},{1,size>1 ? 1 : 0}};
  std::vector<Tensor> v,e;std::vector<int> sites,edgeids;
  tnbp::InitTensorProductState(ctx,p.graph.edges,{2,2},v,sites,owners,e,edgeids,comm);
  auto t=tnbp::MakePauliLayerExponentialTPO<Tensor>(ctx,p,0,0.2);
  std::vector<Tensor> local_v,local_e;
  for(const auto& a:v)local_v.push_back(tcapi::copy(ctx,a));
  for(const auto& a:e)local_e.push_back(tcapi::copy(ctx,a));
  auto local_t=tnbp::MakePauliLayerExponentialTPO<Tensor>(ctx,p,0,.2,sites);
  check(local_t.log_scale==t.log_scale,"local MPI scalar agrees on every rank");
  tnbp::AttachTPO(ctx,p.graph.edges,local_t.tensors,local_v,sites,owners,local_e,edgeids,comm);
  tnbp::AbsorbTPO(ctx,t.tensors,p.graph.edges,v,sites,owners,e,edgeids,comm);
  for(std::size_t i=0;i<v.size();++i)check(tensor_bytes(ctx,v[i])==tensor_bytes(ctx,local_v[i]),"local AttachTPO and global AbsorbTPO give identical V");
  for(std::size_t i=0;i<e.size();++i)check(tensor_bytes(ctx,e[i])==tensor_bytes(ctx,local_e[i]),"local AttachTPO and global AbsorbTPO give identical E");
  // Test-only collection of distributed tensors; production generation has no MPI.
  std::vector<Tensor> all_v(2);
  for(int site=0;site<2;++site) {
    if(rank==owners.at(site)) {
      auto pos=std::find(sites.begin(),sites.end(),site)-sites.begin();
      all_v[site]=tcapi::copy(ctx,v.at(pos));
    }
    tnbp::MpiBcast(ctx,all_v[site],owners.at(site),comm);
  }
  if(rank==0) {
    std::vector<std::vector<C>> flat;
    for(const auto& a:all_v) {
      auto shape=tcapi::shape(ctx,a);std::vector<C> data(tcapi::size(ctx,a));
      tcapi::to_range(ctx,a,data.begin(),[shape](const auto& c){return c[0]+shape[0]*c[1];});flat.push_back(data);
    }
    auto expected=exponential("YZ",0.2*0.37,true);
    for(int bits=0;bits<4;++bits) {
      C value=0;for(int b=0;b<2;++b)value+=flat[0][b+2*(bits&1)]*flat[1][b+2*((bits>>1)&1)];
      check(std::abs(value-expected[bits*4])<1e-13,"AbsorbTPO input/output and Y sign");
    }
  }
}
void absorption_layout(Context& ctx,MPI_Comm comm) {
  int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
  // Unsorted graph, a branching route, an inactive seam and an isolated site.
  auto p=plan({{{.37,0},"IIZYIZ"}},{{1,3},{0,1},{1,2},{3,4}});
  std::map<int,int> owners;
  for(int site:p.graph.sites) owners[site]=size>1 && (site==0 || site==2 || site==4) ? 1 : 0;
  std::vector<Tensor> v,e;std::vector<int> sites,edgeids;
  tnbp::InitTensorProductState(ctx,p.graph.sites,p.graph.edges,{2,2,2,2,2,2},v,sites,owners,e,edgeids,comm);
  std::vector<Tensor> reordered;
  for(auto it=v.rbegin();it!=v.rend();++it)reordered.push_back(tcapi::copy(ctx,*it));
  v.swap(reordered);std::reverse(sites.begin(),sites.end());
  reordered.clear();const auto n=edgeids.size();
  for(std::size_t j=n;j>0;--j)reordered.push_back(tcapi::copy(ctx,e[j-1]));
  for(std::size_t j=n;j>0;--j)reordered.push_back(tcapi::copy(ctx,e[n+j-1]));
  e.swap(reordered);std::reverse(edgeids.begin(),edgeids.end());
  std::vector<Tensor> local_v,local_e;
  for(const auto& a:v)local_v.push_back(tcapi::copy(ctx,a));
  for(const auto& a:e)local_e.push_back(tcapi::copy(ctx,a));
  auto full=tnbp::MakePauliLayerExponentialTPO<Tensor>(ctx,p,0,.2);
  auto local=tnbp::MakePauliLayerExponentialTPO<Tensor>(ctx,p,0,.2,sites);
  check(local.tensors.size()==sites.size() && local.log_scale==full.log_scale,"permuted local layout and remote-root scalar");
  for(int step=1;step<=2;++step) {
    tnbp::AbsorbTPO(ctx,full.tensors,p.graph.sites,p.graph.edges,v,sites,owners,e,edgeids,comm);
    tnbp::AttachTPO(ctx,p.graph.edges,local.tensors,local_v,sites,owners,local_e,edgeids,comm);
    for(std::size_t i=0;i<v.size();++i)check(tensor_bytes(ctx,v[i])==tensor_bytes(ctx,local_v[i]),"permuted Attach/Absorb V byte identity");
    for(std::size_t i=0;i<e.size();++i)check(tensor_bytes(ctx,e[i])==tensor_bytes(ctx,local_e[i]),"permuted Attach/Absorb E byte identity");
    for(std::size_t j=0;j<n;++j) {
      const tnbp::bond_dim_t<Tensor> d=edgeids[j]==3 ? 1 : 1<<step;
      for(auto slot:{j,j+n}) {
        check(tcapi::shape(ctx,local_e[slot])==tnbp::shape_t<Tensor>({d,d}),"active seam expands, unused seam stays unchanged");
        for(tnbp::elem_coor_t<Tensor> a=0;a<d;++a)for(tnbp::elem_coor_t<Tensor> b=0;b<d;++b)
          check(std::abs(tcapi::get_elem(ctx,local_e[slot],{a,b})-C(a==b ? 1/std::sqrt(double(d)):0))<2e-14,"both directed messages have the expected normalized identity");
      }
    }
  }
}
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);int rank;MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  try {
#ifdef TNBP_TEST_CUDA
    int count=0;cudaGetDeviceCount(&count);check(count==1,"one visible GPU per rank");cudaSetDevice(0);
#endif
    Context ctx;tcapi::create_context(ctx);
    planning();operators(ctx);support_weights<Tensor>(ctx);local_generation(ctx);trotter(ctx);single_precision();absorption(ctx,MPI_COMM_WORLD);absorption_layout(ctx,MPI_COMM_WORLD);
    if(rank==0)std::cout<<tnbp::make_timestamp()<<" # pauli_evolution: PASS checks="<<checks
      <<" max_operator_error="<<max_operator_error
      <<" first_order_ratio="<<convergence_ratios[0]
      <<" second_order_ratio="<<convergence_ratios[1]<<'\n';
  } catch(const std::exception& e) {std::cerr<<tnbp::make_timestamp()<<" rank="<<rank<<" "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}
  MPI_Finalize();
}
