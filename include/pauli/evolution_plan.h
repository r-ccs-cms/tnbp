#pragma once
#include "pauli/pauli_string.h"
#include "tnbp/ptns/function/detail/checkpoint_sha256.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace tnbp {
struct PauliHamiltonianTerm {
  double coefficient;
  std::string label;
  std::vector<pauli::LocalOp> local_ops;
};
struct PauliHamiltonian {
  std::size_t num_sites = 0;
  double identity_coefficient = 0;
  std::vector<PauliHamiltonianTerm> terms;
};
struct OperatorGraph {
  std::vector<int> sites;
  std::vector<std::pair<int,int>> edges;
  std::vector<std::vector<std::size_t>> incident_edge_ids;
};
struct PauliRouteSite {
  int site;
  char op;
  std::vector<std::size_t> tree_legs;
};
struct PauliTermRoute {
  std::size_t term_id;
  int root_site;
  std::vector<int> tree_sites;
  std::vector<std::size_t> tree_edge_ids;
  std::vector<PauliRouteSite> local_sites;
};
struct PauliLayer {
  std::vector<std::size_t> term_ids, edge_ids;
};
/** A validated plan is a value object; do not modify its fields after construction. */
struct PauliEvolutionPlan {
  PauliHamiltonian hamiltonian;
  OperatorGraph graph;
  std::vector<PauliTermRoute> routes;
  std::vector<PauliLayer> layers;
  int plan_version = 2;
  std::string plan_digest;
};
struct TrotterAction { std::size_t layer_id; double time_fraction; };

inline PauliHamiltonian NormalizePauliHamiltonian(
    const std::vector<pauli::Term<std::complex<double>>>& input) {
  if(input.empty()) throw std::invalid_argument("Pauli Hamiltonian: empty input cannot define site count");
  PauliHamiltonian h;
  std::map<std::string,std::complex<double>> sums;
  for(const auto& term:input) {
    auto label=term.pauli_string;
    try { pauli::normalize_pauli_label(label); }
    catch(const std::runtime_error& e) { throw std::invalid_argument(e.what()); }
    if(label.empty()) throw std::invalid_argument("Pauli Hamiltonian: empty label");
    if(!h.num_sites) h.num_sites=label.size();
    if(label.size()!=h.num_sites) throw std::invalid_argument("Pauli Hamiltonian: inconsistent label length");
    if(h.num_sites>static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw std::overflow_error("Pauli Hamiltonian: site ID range exceeded");
    if(!std::isfinite(term.coeff.real()) || !std::isfinite(term.coeff.imag()))
      throw std::invalid_argument("Pauli Hamiltonian: nonfinite coefficient for "+label);
    auto& sum=sums[label]; sum+=term.coeff;
    if(!std::isfinite(sum.real()) || !std::isfinite(sum.imag()))
      throw std::overflow_error("Pauli Hamiltonian: coefficient sum overflow for "+label);
  }
  for(const auto& [label,c]:sums) {
    if(c.imag()!=0) throw std::invalid_argument("Pauli Hamiltonian: nonreal coefficient for "+label);
    if(c.real()==0) continue;
    auto ops=pauli::string_to_non_identity_localops(label);
    if(ops.empty()) h.identity_coefficient=c.real();
    else h.terms.push_back({c.real(),label,std::move(ops)});
  }
  return h;
}
namespace pauli_plan_detail {
inline bool commute(const std::string& a,const std::string& b) {
  bool odd=false;
  for(std::size_t i=0;i<a.size();++i)
    if(a[i]!='I' && b[i]!='I' && a[i]!=b[i]) odd=!odd;
  return !odd;
}
// Shortest paths from one support site. Sorted adjacency makes the first path
// to each vertex lexicographically first among shortest paths to that vertex.
inline std::vector<int> path_to_tree(int src,const std::set<int>& tree,
    const std::vector<std::vector<int>>& adjacency) {
  std::vector<int> previous(adjacency.size(),-1), distance(adjacency.size(),-1);
  std::queue<int> pending; pending.push(src); distance[src]=0;
  std::vector<int> best;
  while(!pending.empty()) {
    int v=pending.front(); pending.pop();
    if(tree.count(v)) {
      std::vector<int> path;
      for(int u=v;u!=-1;u=previous[u]) path.push_back(u);
      std::reverse(path.begin(),path.end());
      if(best.empty() || std::make_tuple(path.size(),v,path)<std::make_tuple(best.size(),best.back(),best))
        best=std::move(path);
      continue; // Paths stop at their first tree vertex.
    }
    for(int u:adjacency[v]) if(distance[u]<0) {
      distance[u]=distance[v]+1; previous[u]=v; pending.push(u);
    }
  }
  return best;
}
inline std::string digest(const PauliEvolutionPlan& p) {
  std::ostringstream s; s.imbue(std::locale::classic());
  s<<"TNBP_PAULI_PLAN "<<p.plan_version<<'\n'<<p.hamiltonian.num_sites<<' '
   <<std::hexfloat<<p.hamiltonian.identity_coefficient<<'\n';
  for(const auto& t:p.hamiltonian.terms) s<<t.label<<' '<<t.coefficient<<'\n';
  s<<"EDGES\n";
  for(auto e:p.graph.edges) s<<e.first<<' '<<e.second<<'\n';
  s<<"ROUTES\n";
  for(const auto& r:p.routes) {
    s<<r.term_id<<' '<<r.root_site<<":";
    for(auto e:r.tree_edge_ids) s<<' '<<e;
    s<<'\n';
  }
  s<<"LAYERS\n";
  for(const auto& l:p.layers) {for(auto t:l.term_ids) s<<t<<' '; s<<'\n';}
  return portable_detail::sha256(s.str());
}
} // namespace pauli_plan_detail

inline PauliEvolutionPlan MakePauliEvolutionPlan(const PauliHamiltonian& h,
    const std::vector<std::pair<int,int>>& ordered_edges) {
  // Revalidate public input rather than relying on callers to normalize it.
  if(!h.num_sites || h.num_sites>static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("Pauli plan: invalid site count");
  std::vector<pauli::Term<std::complex<double>>> raw;
  raw.push_back({{h.identity_coefficient,0},std::string(h.num_sites,'I')});
  for(const auto& t:h.terms) {
    if(t.label.size()!=h.num_sites) throw std::invalid_argument("Pauli plan: inconsistent term label");
    raw.push_back({{t.coefficient,0},t.label});
  }
  PauliEvolutionPlan p; p.hamiltonian=NormalizePauliHamiltonian(raw);
  p.graph.edges=ordered_edges; p.graph.sites.resize(h.num_sites);
  std::iota(p.graph.sites.begin(),p.graph.sites.end(),0);
  p.graph.incident_edge_ids.resize(h.num_sites);
  std::vector<std::vector<int>> adjacency(h.num_sites);
  std::map<std::pair<int,int>,std::size_t> edge_ids;
  for(std::size_t e=0;e<ordered_edges.size();++e) {
    auto [a,b]=ordered_edges[e];
    if(a<0 || b<0 || static_cast<std::size_t>(a)>=h.num_sites || static_cast<std::size_t>(b)>=h.num_sites || a==b)
      throw std::invalid_argument("Pauli plan: invalid edge "+std::to_string(e));
    if(!edge_ids.emplace(std::minmax(a,b),e).second)
      throw std::invalid_argument("Pauli plan: duplicate undirected edge "+std::to_string(e));
    adjacency[a].push_back(b); adjacency[b].push_back(a);
    p.graph.incident_edge_ids[a].push_back(e); p.graph.incident_edge_ids[b].push_back(e);
  }
  for(auto& a:adjacency) std::sort(a.begin(),a.end());
  for(std::size_t id=0;id<p.hamiltonian.terms.size();++id) {
    const auto& term=p.hamiltonian.terms[id];
    int root=static_cast<int>(term.local_ops.front().qubit);
    std::set<int> tree{root}; std::set<std::size_t> used;
    while(true) {
      std::vector<int> best;
      for(const auto& op:term.local_ops) if(!tree.count(static_cast<int>(op.qubit))) {
        auto path=pauli_plan_detail::path_to_tree(static_cast<int>(op.qubit),tree,adjacency);
        if(path.empty()) throw std::invalid_argument("Pauli plan: disconnected support for term "+term.label);
        if(best.empty() || std::make_tuple(path.size(),path.front(),path.back(),path)<
            std::make_tuple(best.size(),best.front(),best.back(),best)) best=std::move(path);
      }
      if(best.empty()) break;
      for(std::size_t k=1;k<best.size();++k) used.insert(edge_ids.at(std::minmax(best[k-1],best[k])));
      tree.insert(best.begin(),best.end());
    }
    if(used.size()+1!=tree.size()) throw std::logic_error("Pauli plan: route is not a tree");
    PauliTermRoute route{id,root,{tree.begin(),tree.end()},{used.begin(),used.end()},{}};
    for(int v:tree) {
      PauliRouteSite local{v,term.label[h.num_sites-1-static_cast<std::size_t>(v)],{}};
      const auto& incident=p.graph.incident_edge_ids[v];
      for(std::size_t k=0;k<incident.size();++k) if(used.count(incident[k])) local.tree_legs.push_back(k);
      route.local_sites.push_back(std::move(local));
    }
    p.routes.push_back(std::move(route));
    std::size_t layer=0;
    for(;layer<p.layers.size();++layer) {
      bool fits=true;
      for(auto other:p.layers[layer].term_ids)
        if(!pauli_plan_detail::commute(term.label,p.hamiltonian.terms[other].label)) {fits=false;break;}
      // Exclude every route vertex, including intermediate identity sites.
      if(fits) for(auto other:p.layers[layer].term_ids) {
        const auto& occupied=p.routes[other].tree_sites;
        if(std::any_of(occupied.begin(),occupied.end(),[&](int site){return tree.count(site)!=0;})) {
          fits=false;break;
        }
      }
      if(fits) break;
    }
    if(layer==p.layers.size()) p.layers.push_back({{}, {}});
    auto& l=p.layers[layer]; l.term_ids.push_back(id);
    l.edge_ids.insert(l.edge_ids.end(),used.begin(),used.end()); std::sort(l.edge_ids.begin(),l.edge_ids.end());
  }
  p.plan_digest=pauli_plan_detail::digest(p);
  return p;
}
/** One edge per term tree per batch, children before parent; no MPI or tensors.
 * Ready edges within a tree use global edge ID order. Trees are visited in
 * term ID order; a conflicting head waits, never skipping ahead in that tree.
 */
inline std::vector<std::vector<int>> MakePauliLayerTruncationBatches(
    const PauliEvolutionPlan& p,std::size_t layer) {
  if(layer>=p.layers.size()) throw std::invalid_argument("Pauli batches: layer ID out of range");
  if(p.graph.edges.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::overflow_error("Pauli batches: edge ID range exceeded");
  auto terms=p.layers[layer].term_ids;std::sort(terms.begin(),terms.end());
  std::vector<std::vector<int>> orders;
  for(auto term:terms) {
    const auto& route=p.routes.at(term);
    std::map<int,std::vector<std::pair<int,int>>> adjacency;
    for(auto id:route.tree_edge_ids) {
      auto [a,b]=p.graph.edges.at(id);
      adjacency[a].push_back({b,static_cast<int>(id)});
      adjacency[b].push_back({a,static_cast<int>(id)});
    }
    std::map<int,int> parent,parent_edge,children;
    parent[route.root_site]=route.root_site;
    std::queue<int> pending;pending.push(route.root_site);
    while(!pending.empty()) {
      int v=pending.front();pending.pop();
      for(auto [u,e]:adjacency[v]) if(!parent.count(u)) {
        parent[u]=v;parent_edge[u]=e;++children[v];pending.push(u);
      }
    }
    std::set<std::pair<int,int>> ready; // edge ID, child site
    for(auto [v,e]:parent_edge) if(children[v]==0) ready.insert({e,v});
    std::vector<int> order;
    while(!ready.empty()) {
      auto [e,v]=*ready.begin();ready.erase(ready.begin());order.push_back(e);
      int up=parent.at(v);
      if(--children[up]==0 && up!=route.root_site) ready.insert({parent_edge.at(up),up});
    }
    if(order.size()!=route.tree_edge_ids.size()) throw std::logic_error("Pauli batches: invalid route tree");
    orders.push_back(std::move(order));
  }
  std::vector<std::size_t> next(orders.size());std::vector<std::vector<int>> batches;
  for(;;) {
    std::set<int> occupied;std::vector<int> batch;
    for(std::size_t t=0;t<orders.size();++t) if(next[t]<orders[t].size()) {
      int e=orders[t][next[t]];auto [a,b]=p.graph.edges[e];
      if(occupied.count(a) || occupied.count(b)) continue;
      occupied.insert(a);occupied.insert(b);batch.push_back(e);++next[t];
    }
    if(batch.empty()) break;
    std::sort(batch.begin(),batch.end());batches.push_back(std::move(batch));
  }
  return batches;
}

inline std::vector<TrotterAction> MakeTrotterActions(const PauliEvolutionPlan& p,int order=2) {
  if(order!=1 && order!=2) throw std::invalid_argument("Pauli Trotter: order must be 1 or 2");
  std::vector<TrotterAction> actions;
  auto n=p.layers.size();
  if(!n) return actions;
  if(order==1) {for(std::size_t i=0;i<n;++i) actions.push_back({i,1});}
  else {
    for(std::size_t i=0;i+1<n;++i) actions.push_back({i,0.5});
    actions.push_back({n-1,1});
    for(std::size_t i=n-1;i>0;--i) actions.push_back({i-1,0.5});
  }
  return actions;
}
} // namespace tnbp
