#pragma once
#include "pauli/evolution_plan.h"
#include "tnbp/framework/typedef.h"
#include <array>
#include <type_traits>

namespace tnbp {
// Numerical layout identifier used by callers that persist evolution settings.
inline constexpr const char* PauliExponentialWeights = "pauli-support-equal-v1";
template<class TenT> struct ScaledPauliTPO {
  std::vector<TenT> tensors;
  double log_scale = 0; // Exact operator = exp(log_scale) * represented TPO.
};
namespace pauli_tpo_detail {
using Matrix=std::array<std::complex<double>,4>; // row=output, column=input
inline Matrix matrix(char op) {
  using C=std::complex<double>;
  switch(op) {
    case 'I': return {1,0,0,1};
    case 'X': return {0,1,1,0};
    case 'Y': return {0,C(0,-1),C(0,1),0};
    case 'Z': return {1,0,0,-1};
  }
  throw std::invalid_argument("Pauli TPO: invalid local operator");
}
inline Matrix multiply(const Matrix& a,const Matrix& b) {
  Matrix c{};
  for(int o=0;o<2;++o) for(int i=0;i<2;++i) for(int k=0;k<2;++k) c[2*o+i]+=a[2*o+k]*b[2*k+i];
  return c;
}
inline std::pair<double,double> coefficients(double x) {
  // Avoid overflow even in forming -2*abs(x), and cancellation near zero.
  const double difference=std::abs(x)>std::numeric_limits<double>::max()/2 ? 1 : -std::expm1(-2*std::abs(x));
  return {1-0.5*difference,-std::copysign(0.5*difference,x)};
}

template<class TenT>
ScaledPauliTPO<TenT> generate(context_handle_t<TenT>& ctx,const PauliEvolutionPlan& p,
    const std::vector<std::size_t>& ids,double dt,const std::vector<int>& target_sites) {
  using Elem=elem_t<TenT>;
  static_assert(std::is_same_v<Elem,std::complex<float>> || std::is_same_v<Elem,std::complex<double>>,
                "Pauli exponential TPO requires complex float/double tensors");
  if(!std::isfinite(dt) || dt<0) throw std::invalid_argument("Pauli TPO: time step must be finite and nonnegative");
  std::set<int> seen;
  for(int site:target_sites)
    if(site<0 || static_cast<std::size_t>(site)>=p.graph.sites.size() || !seen.insert(site).second)
      throw std::invalid_argument("Pauli TPO: invalid or duplicate target site");
  ScaledPauliTPO<TenT> result;
  result.tensors.reserve(target_sites.size());
  std::vector<std::pair<double,double>> weights;
  std::vector<std::pair<double,double>> site_weights;
  std::vector<bool> active(p.graph.edges.size(),false);
  for(auto id:ids) {
    const double x=dt*p.hamiltonian.terms.at(id).coefficient;
    if(!std::isfinite(x)) throw std::overflow_error("Pauli TPO: coefficient times time step overflow");
    result.log_scale+=std::abs(x);
    if(!std::isfinite(result.log_scale)) throw std::overflow_error("Pauli TPO: log scale overflow");
    const auto weight=coefficients(x);
    weights.push_back(weight);
    const auto support=p.hamiltonian.terms[id].local_ops.size();
    const double power=1.0/static_cast<double>(support);
    site_weights.emplace_back(std::pow(weight.first,power),
                              std::pow(std::abs(weight.second),power));
    if(dt!=0) for(auto e:p.routes.at(id).tree_edge_ids) {
      if(active[e]) throw std::invalid_argument("Pauli TPO: repeated layer edge");
      active[e]=true;
    }
  }
  for(int v:target_sites) {
    const auto& incident=p.graph.incident_edge_ids[v];
    if(incident.size()>static_cast<std::size_t>(std::numeric_limits<order_t<TenT>>::max())-2)
      throw std::overflow_error("Pauli TPO: tensor order exceeds backend range");
    shape_t<TenT> shape(incident.size()+2,1);
    for(std::size_t k=0;k<incident.size();++k) shape[k]=active[incident[k]] ? 2 : 1;
    shape[incident.size()]=2; shape[incident.size()+1]=2;
    std::size_t count=1;
    for(auto d:shape) {
      if(count>std::numeric_limits<std::size_t>::max()/static_cast<std::size_t>(d))
        throw std::overflow_error("Pauli TPO: tensor size overflow at site "+std::to_string(v));
      count*=static_cast<std::size_t>(d);
    }
    std::vector<Elem> data(count);
    for(std::size_t linear=0;linear<count;++linear) {
      elem_coors_t<TenT> coord(shape.size()); auto remainder=linear;
      for(std::size_t k=0;k<shape.size();++k) {coord[k]=remainder%shape[k]; remainder/=shape[k];}
      auto product=matrix('I'); bool nonzero=true;
      if(dt!=0) for(std::size_t j=0;j<ids.size();++j) {
        const auto& route=p.routes[ids[j]];
        auto local=std::lower_bound(route.local_sites.begin(),route.local_sites.end(),v,
            [](const PauliRouteSite& s,int site){return s.site<site;});
        if(local==route.local_sites.end() || local->site!=v) continue;
        Matrix factor{};
        const auto [a,b]=weights[j];
        if(local->tree_legs.empty()) { // Single-site term: no virtual branch.
          const auto identity=matrix('I'), op=matrix(local->op);
          for(int k=0;k<4;++k) factor[k]=a*identity[k]+b*op[k];
        } else {
          auto branch=coord[local->tree_legs.front()];
          for(auto leg:local->tree_legs) if(coord[leg]!=branch) nonzero=false;
          if(!nonzero) break;
          factor=matrix(branch ? local->op : 'I');
          // Share each branch coefficient once per non-identity Pauli site.
          // Routing I sites, including branch points, remain pure copy tensors.
          if(local->op!='I') {
            const auto [local_a,local_b]=site_weights[j];
            double weight=branch ? local_b : local_a;
            // Keep a negative coefficient's sign at one deterministic site.
            if(branch && b<0 && v==route.root_site) weight=-weight;
            for(auto& value:factor) value*=weight;
          }
        }
        product=multiply(factor,product);
      }
      if(nonzero) data[linear]=static_cast<Elem>(product[2*coord[incident.size()+1]+coord[incident.size()]]);
    }
    result.tensors.push_back(tcapi::assign_from_range<TenT>(ctx,shape,data.begin(),[&shape](const auto& c) {
      std::size_t address=0,stride=1;
      for(std::size_t k=0;k<shape.size();++k) {address+=static_cast<std::size_t>(c[k])*stride; stride*=shape[k];}
      return address;
    }));
  }
  return result;
}
} // namespace pauli_tpo_detail
/** Construct one scaled exponential on target_sites, preserving their order.
 * log_scale describes the complete operator even for an empty site list.
 * No MPI or logging; target_sites can be a rank's local site_idx.
 */
template<class TenT>
ScaledPauliTPO<TenT> MakePauliTermExponentialTPO(context_handle_t<TenT>& ctx,
    const PauliEvolutionPlan& p,std::size_t id,double dt,const std::vector<int>& target_sites) {
  if(id>=p.hamiltonian.terms.size()) throw std::invalid_argument("Pauli TPO: term ID out of range");
  return pauli_tpo_detail::generate<TenT>(ctx,p,{id},dt,target_sites);
}
/** Compatibility overload: construct all sites in plan.graph.sites order. */
template<class TenT>
ScaledPauliTPO<TenT> MakePauliTermExponentialTPO(context_handle_t<TenT>& ctx,
    const PauliEvolutionPlan& p,std::size_t id,double dt) {
  return MakePauliTermExponentialTPO<TenT>(ctx,p,id,dt,p.graph.sites);
}
/** Construct a planned layer on target_sites. The identity Hamiltonian term
 * is excluded. log_scale is global and must not be summed over MPI ranks.
 */
template<class TenT>
ScaledPauliTPO<TenT> MakePauliLayerExponentialTPO(context_handle_t<TenT>& ctx,
    const PauliEvolutionPlan& p,std::size_t layer,double dt,const std::vector<int>& target_sites) {
  if(layer>=p.layers.size()) throw std::invalid_argument("Pauli TPO: layer ID out of range");
  return pauli_tpo_detail::generate<TenT>(ctx,p,p.layers[layer].term_ids,dt,target_sites);
}
/** Compatibility overload: construct all sites in plan.graph.sites order. */
template<class TenT>
ScaledPauliTPO<TenT> MakePauliLayerExponentialTPO(context_handle_t<TenT>& ctx,
    const PauliEvolutionPlan& p,std::size_t layer,double dt) {
  return MakePauliLayerExponentialTPO<TenT>(ctx,p,layer,dt,p.graph.sites);
}
} // namespace tnbp
