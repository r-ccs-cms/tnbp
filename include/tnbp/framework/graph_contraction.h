/// This file is a part of r-ccs-cms/tnbp
/**
 * @file tnbp/framework/graph_contraction.h
 * @brief Graph-aware tensor contraction and compatible graph updates.
 */
#ifndef TNBP_FRAMEWORK_GRAPH_CONTRACTION_H
#define TNBP_FRAMEWORK_GRAPH_CONTRACTION_H

#include "tnbp/framework/graph_contraction_plan.h"
#include "tnbp/framework/mpiutility.h"
#include "tnbp/framework/typedef.h"
#include <cstddef>
#include <iterator>
#include <map>
#include <set>
#include <limits>
#include <type_traits>
#include <stdexcept>
#include <numeric>
#include <utility>

#include <vector>
#include <algorithm>

namespace tnbp {

namespace detail {
// Partial environment reform explicitly marks unavailable graph sites with -1.
// Both operands must either have tensors or both be graph-only. Missing shared
// neighbors are permitted only in this internal mode; their tensors are reformed
// separately by the sampling caller before they are combined.
enum class GraphContractionMode { complete, partial_environment };

// The temporary factor graph used by both tensor and graph-only bond shifts.
template<class IntT>
std::vector<std::pair<IntT,IntT>> make_graph_shift_edges(
    const std::vector<std::pair<IntT,IntT>>& edges, IntT source, IntT destination,
    const std::set<IntT>& selected, IntT temporary) {
    auto updated=edges;
    for(auto& edge:updated) {
        if(edge.first!=source && edge.second!=source) continue;
        const auto other=edge.first==source ? edge.second : edge.first;
        if(other==destination || selected.count(other)) edge=make_edge(temporary,other);
    }
    updated.push_back(make_edge(temporary,source));
    if(find_edge(edges,source,destination)==edges.end())
        updated.push_back(make_edge(temporary,destination));
    return updated;
}

// Collective preflight for graph contraction, bond shifts and BMPS sweeps.
// Every tensor is indexed by current incident-edge order. Only required sites
// need tensors; graph-only external sites may be absent from owners/site_idx.
template<class TenT, class IntT>
std::map<IntT,std::size_t> validate_graph_tensors(
    context_handle_t<TenT>& ctx, const std::vector<IntT>& required,
    const std::vector<std::pair<IntT,IntT>>& edges,
    const std::vector<TenT>& tensors,const std::vector<IntT>& site_idx,
    const std::map<IntT,IntT>& owners,MPI_Comm comm, bool invalid_operation=false) {
    int rank,size; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&size);
    std::map<IntT,std::size_t> required_sites;
    std::set<std::pair<IntT,IntT>> unique_edges;
    int invalid=invalid_operation ? 1 : 0;
    for(std::size_t j=0;j<required.size();++j) {
      if(!required_sites.emplace(required[j],j).second) invalid=1;
      const auto owner=owners.find(required[j]);
      if(owner==owners.end() || owner->second<0 || owner->second>=size) invalid=1;
    }
    for(const auto& [a,b]:edges) {
      if(a==b || !unique_edges.insert(make_edge(a,b)).second) invalid=1;
    }
    if(tensors.size()!=site_idx.size() || edges.size()>static_cast<std::size_t>(std::numeric_limits<int>::max())) invalid=1;
    std::map<IntT,std::size_t> local;
    for(std::size_t j=0;j<site_idx.size();++j) {
      if(!local.emplace(site_idx[j],j).second) invalid=1;
      const auto owner=owners.find(site_idx[j]);
      if(owner==owners.end() || owner->second!=rank) invalid=1;
    }
    // Compare dimensions at both endpoints of every internal required edge, also
    // when those endpoints live on different ranks.
    std::vector<unsigned long long> dim_min(edges.size(),std::numeric_limits<unsigned long long>::max());
    std::vector<unsigned long long> dim_max(edges.size(),0);
    try {
      for(auto site:required) {
        const auto owner=owners.find(site);
        if(owner==owners.end()) continue;
        const auto it=local.find(site);
        if((it!=local.end()) != (owner->second==rank)) {invalid=1;continue;}
        if(it==local.end() || it->second>=tensors.size()) continue;
        const auto bonds=GetSurroundingBondIndex(site,edges);
        const auto shape=tcapi::shape(ctx,tensors[it->second]);
        if(shape.size()!=bonds.size()) {invalid=1;continue;}
        for(std::size_t k=0;k<shape.size();++k) {
          if(shape[k]<1) {invalid=1;continue;}
          const auto e=static_cast<std::size_t>(bonds[k]);
          const auto d=static_cast<unsigned long long>(shape[k]);
          dim_min[e]=std::min(dim_min[e],d);dim_max[e]=std::max(dim_max[e],d);
        }
      }
    } catch(const std::exception&) {invalid=1;}
    int any_invalid=0;
    MPI_Allreduce(&invalid,&any_invalid,1,MPI_INT,MPI_MAX,comm);
    if(any_invalid) throw std::invalid_argument("graph operation: invalid graph, required sites, owner or tensor axes");
    if(!edges.empty()) {
      MPI_Allreduce(MPI_IN_PLACE,dim_min.data(),static_cast<int>(edges.size()),MPI_UNSIGNED_LONG_LONG,MPI_MIN,comm);
      MPI_Allreduce(MPI_IN_PLACE,dim_max.data(),static_cast<int>(edges.size()),MPI_UNSIGNED_LONG_LONG,MPI_MAX,comm);
      for(std::size_t e=0;e<edges.size();++e)
        if(required_sites.count(edges[e].first) && required_sites.count(edges[e].second) && dim_min[e]!=dim_max[e])
          throw std::invalid_argument("graph operation: internal bond dimensions disagree");
    }
    return local;
}

// Fuse removed_axis into kept_axis in that order. Other axes retain their
// relative order; the fused axis occupies kept_axis's position after removal.
// Both ends of a fused graph bond must call this with the same axis roles.
template<class TenT>
void fuse_graph_bond_axes(context_handle_t<TenT>& ctx, TenT& tensor,
                          std::size_t kept_axis, std::size_t removed_axis) {
    const auto shape=tcapi::shape(ctx,tensor);
    if(kept_axis==removed_axis || kept_axis>=shape.size() || removed_axis>=shape.size())
        throw std::invalid_argument("graph bond fusion: invalid axes");
    tcapi::List<bond_idx_t<TenT>> permutation;
    shape_t<TenT> fused_shape;
    for(std::size_t axis=0;axis<shape.size();++axis) {
        if(axis==removed_axis) continue;
        permutation.push_back(static_cast<bond_idx_t<TenT>>(axis));
        if(axis==kept_axis) {
            permutation.push_back(static_cast<bond_idx_t<TenT>>(removed_axis));
            fused_shape.push_back(shape[kept_axis]*shape[removed_axis]);
        } else {
            fused_shape.push_back(shape[axis]);
        }
    }
    tcapi::transpose(ctx,tensor,permutation);
    tcapi::reshape(ctx,tensor,fused_shape);
}
} // namespace detail

  template<class IntT>
  void graph_rename_site(std::vector<std::pair<IntT,IntT>>& edges,
                         std::vector<IntT>& ids,std::map<IntT,IntT>& owners,
                         IntT source,IntT destination) {
    if(!owners.count(source)) throw std::invalid_argument("rename source is absent");
    if(source==destination) return;
    if(owners.count(destination)) throw std::invalid_argument("rename destination already exists");
    for(auto& edge:edges) {
      if(edge.first==source) edge.first=destination;
      if(edge.second==source) edge.second=destination;
      edge=make_edge(edge.first,edge.second);
    }
    for(auto& id:ids) if(id==source) id=destination;
    auto node=owners.extract(source);node.key()=destination;owners.insert(std::move(node));
  }

  /**
   * @brief Compute the symmetric difference of two label sets.
   *
   * Given two label lists `label_a` and `label_b`, this function computes
   * the set of elements that appear in either list but not in both
   * (i.e., the symmetric difference), and stores the result in `label_c`.
   *
   * The output `label_c` is sorted in ascending order.
   *
   * @tparam BondLabelT Integer type of labels.
   *
   * @param[in]  label_a First label list (no duplicates).
   * @param[in]  label_b Second label list (no duplicates).
   * @param[out] label_c Output label list containing the symmetric difference,
   *                     sorted in ascending order.
   *
   * @pre
   * - `label_a` and `label_b` contain no duplicate elements.
   *
   * @note
   * - The input lists are not modified.
   * - Internally, copies of the inputs are sorted before computing the result.
   */
  template <typename BondLabelT>
  void tensor_contraction_label_helper(const tcapi::List<BondLabelT> & label_a,
				       const tcapi::List<BondLabelT> & label_b,
				       tcapi::List<BondLabelT> & label_c) {
    tcapi::List<BondLabelT> a = label_a;
    tcapi::List<BondLabelT> b = label_b;
    std::sort(a.begin(),a.end());
    std::sort(b.begin(),b.end());

    label_c.clear();

    std::set_symmetric_difference(
	 a.begin(),a.end(),
	 b.begin(),b.end(),
	 std::back_inserter(label_c));
  }


namespace detail {
  template <typename TenT, typename IntT>
  void graph_tensor_contraction_impl(context_handle_t<TenT> & ctx,
				std::vector<std::pair<IntT,IntT>> & edges,
				std::vector<TenT> & W,
				std::vector<IntT> & site_idx,
				std::map<IntT,IntT> & site_to_mpi_rank,
				IntT site_a, IntT site_b, IntT site_c,
				MPI_Comm comm, GraphContractionMode mode) {

    using BondLabelT = typename tcapi::tensor_traits<TenT>::bond_label_t;

    // All ranks validate before any point-to-point transfer or tensor mutation.
    // Unaffected external sites may remain absent from tensors and owners.
    GraphContractionPlan<IntT> plan;
    bool invalid=site_c!=site_a && site_to_mpi_rank.count(site_c);
    try {plan=make_graph_contraction_plan(edges,site_a,site_b,site_c);}
    catch(const std::invalid_argument&) {invalid=true;}
    std::vector<IntT> candidates{site_a,site_b};
    for(const auto& fusion:plan.fusions) candidates.push_back(fusion.neighbor);
    std::vector<IntT> required;
    const auto graph_only=[&](IntT site) {
      const auto owner=site_to_mpi_rank.find(site);
      return mode==GraphContractionMode::partial_environment &&
             owner!=site_to_mpi_rank.end() && owner->second==IntT(-1);
    };
    // A graph-only reform never sends a real tensor to a missing owner.
    invalid=invalid || graph_only(site_a)!=graph_only(site_b);
    for(auto site:candidates) if(!graph_only(site)) required.push_back(site);
    validate_graph_tensors(ctx,required,edges,W,site_idx,site_to_mpi_rank,comm,invalid);

    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);

    tcapi::List<BondLabelT> label_a(plan.bonds_a.begin(),plan.bonds_a.end());
    tcapi::List<BondLabelT> label_b(plan.bonds_b.begin(),plan.bonds_b.end());
    tcapi::List<BondLabelT> label_c(plan.result_bonds.begin(),plan.result_bonds.end());

    int mpi_rank_a = static_cast<int>(site_to_mpi_rank.at(site_a));
    int mpi_rank_b = static_cast<int>(site_to_mpi_rank.at(site_b));
    int mpi_type = 0;
    if( mpi_rank == mpi_rank_a ) {
      mpi_type += 2;
    }
    if( mpi_rank == mpi_rank_b ) {
      mpi_type += 1;
    }

    if( mpi_type > 0 ) {
      TenT B;
      if( mpi_type == 1 || mpi_type == 3 ) {
	auto it_adrs_b = std::find(site_idx.begin(),site_idx.end(),site_b);
	auto adrs_b = std::distance(site_idx.begin(),it_adrs_b);
	B = tcapi::copy(ctx,W[adrs_b]);
      }
      if( mpi_type == 1 ) {
	MpiSend(ctx,B,mpi_rank_a,comm);
      }
      if( mpi_type == 2 ) {
	MpiRecv(ctx,B,mpi_rank_b,comm);
      }
      if( mpi_type == 2 || mpi_type == 3 ) {
	auto it_adrs_a = std::find(site_idx.begin(),site_idx.end(),site_a);
	auto adrs_a = std::distance(site_idx.begin(),it_adrs_a);
#ifdef HELPER_DEBUG
	auto shape_a = tcapi::shape(ctx,W[adrs_a]);
	std::cout << " graph_tensor_contraction: Tensor A at site " << site_a << ": shape =";
	for(const auto & dim : shape_a) {
	  std::cout << " " << dim;
	}
	std::cout << ", label = ";
	for(const auto & lbl : label_a) {
	  std::cout << " " << lbl;
	}
	std::cout << std::endl;
	auto shape_b = tcapi::shape(ctx,B);
	std::cout << " graph_tensor_contraction: Tensor B at site " << site_b << ": shape =";
	for(const auto & dim : shape_b) {
	  std::cout << " " << dim;
	}
	std::cout << ", label = ";
	for(const auto & lbl : label_b) {
	  std::cout << " " << lbl;
	}
	std::cout << std::endl;
#endif
	tcapi::contract(ctx,W[adrs_a],label_a,B,label_b,W[adrs_a],label_c);
#ifdef HELPER_DEBUG
	auto shape_c = tcapi::shape(ctx,W[adrs_a]);
	std::cout << " graph_tensor_contraction: Tensor C labeled as " << site_c << ": shape =";
	for(const auto & dim : shape_c) {
	  std::cout << " " << dim;
	}
	std::cout << ", label = ";
	for(const auto & lbl : label_c) {
	  std::cout << " " << lbl;
	}
	std::cout << std::endl;
#endif
      }
    }

    auto result_bonds=plan.result_bonds;
    for(const auto& fusion:plan.fusions) {
      if(mpi_rank==mpi_rank_a) {
        const auto tensor=std::find(site_idx.begin(),site_idx.end(),site_a)-site_idx.begin();
        const auto kept=std::find(result_bonds.begin(),result_bonds.end(),fusion.kept_edge);
        const auto removed=std::find(result_bonds.begin(),result_bonds.end(),fusion.removed_edge);
        fuse_graph_bond_axes(ctx,W[tensor],kept-result_bonds.begin(),removed-result_bonds.begin());
        result_bonds.erase(removed);
      }
      if(mpi_rank==site_to_mpi_rank.at(fusion.neighbor)) {
        const auto tensor=std::find(site_idx.begin(),site_idx.end(),fusion.neighbor)-site_idx.begin();
        const auto bonds=GetSurroundingBondIndex(fusion.neighbor,edges);
        fuse_graph_bond_axes(ctx,W[tensor],
          std::find(bonds.begin(),bonds.end(),fusion.kept_edge)-bonds.begin(),
          std::find(bonds.begin(),bonds.end(),fusion.removed_edge)-bonds.begin());
      }
    }
    edges=std::move(plan.edges);

    auto it_site_b = std::find(site_idx.begin(),site_idx.end(),site_b);
    if( it_site_b != site_idx.end() ) {
      auto adrs_b = std::distance(site_idx.begin(),it_site_b);
      auto it_w_b = W.begin();
      std::advance(it_w_b,adrs_b);
      W.erase(it_w_b);
      site_idx.erase(it_site_b);
    }
    auto it_site_a = std::find(site_idx.begin(),site_idx.end(),site_a);
    if( it_site_a != site_idx.end() ) {
      *it_site_a = site_c;
    }

    auto node = site_to_mpi_rank.extract(site_a);
    if (!node.empty()) {
      node.key() = site_c;
      site_to_mpi_rank.insert(std::move(node));
    }

    // Remove ownership of the contracted-away site.
    auto it_site_b_rank = site_to_mpi_rank.find(site_b);
    if( it_site_b_rank != site_to_mpi_rank.end() ) {
      site_to_mpi_rank.erase(it_site_b_rank);
    }
  }

} // namespace detail

  /**
   * @brief A function for performing tensor contraction between neighboring nodes on a graph in a
   *        distributed environment, with automatic updates to the graph structure and associated tensors.
   * @pre edges describes the current graph, but W/site_idx need not contain
   *      every site in it. On each owning MPI rank they must contain site_a,
   *      site_b, and every common neighbor of site_a and site_b. Common-neighbor
   *      tensors are transposed/reshaped when the two incident bonds are fused.
   *      Other sites (including neighbors of only one operand) need no tensor.
   *      Required sites must have valid entries in site_to_mpi_rank, and tensor
   *      axes must match the incident-edge order in edges.
   * @pre site_a and site_b are distinct adjacent sites. site_c must be site_a
   *      or a fresh label, not site_b or another existing site.
   * @note Collective on comm: all ranks must agree on edges, arguments and owners.
   *       Required tensors, axes, internal bond dimensions and structural inputs
   *       are checked collectively before mutation. Backend failures after this
   *       preflight may still leave partial output.
   * @note Ownership is asymmetric: the result stays on the original owner of
   *       site_a, regardless of the label site_c. site_b is transferred there
   *       when needed. Swapping site_a/site_b changes the destination rank.
   *       Use site_c == site_a or a fresh label: site_c == site_b is not supported
   *       by the current ownership-map update.
   * @tparam TenT type of tensor
   * @tparam IntT type of site label
   *
   * @param[in/out] ctx Context handler for tensor computing interface.
   * @param[in/out] edges Edges to define the graph. On the output, edges is updated after taking contraction.
   * @param[in/out] W Tensor W[i] corresponds to the tensor at site labeled by site_idx[i].
   * @param[in/out] site_idx tensors site label at each mpi rank. The size of site_idx is identical to the size of W.
   * @param[in/out] site_to_mpi_rank A global map to find which site tensor is handled by which node.
   * @param[in] site_a W[i] with site_idx[i]=site_a will be contracted, and site_idx[i] will be replaced as site_c.
   * @param[in] site_b W[i] with site_idx[i]=site_b will be contracted.
   * @param[in] site_c Result label; inherits the original owner of site_a.
   *                  site_idx and site_to_mpi_rank are updated accordingly.
   *
   */
// Public contraction requires both operands and all affected common neighbors.
template<class TenT, class IntT>
void graph_tensor_contraction(context_handle_t<TenT>& ctx,
    std::vector<std::pair<IntT,IntT>>& edges, std::vector<TenT>& tensors,
    std::vector<IntT>& site_idx, std::map<IntT,IntT>& owners,
    IntT site_a, IntT site_b, IntT site_c, MPI_Comm comm) {
    detail::graph_tensor_contraction_impl(ctx,edges,tensors,site_idx,owners,
        site_a,site_b,site_c,comm,detail::GraphContractionMode::complete);
}


/**
 * @brief Move selected source--other bonds onto destination by a full SVD.
 *
 * For each other in bond_neighbors, source--other becomes destination--other.
 * A source--destination SVD bond always remains. If that bond already existed,
 * it is included in the moved factor and replaced by the new SVD bond. Otherwise
 * a temporary dimension-one bond implements the outer product at destination.
 * Existing destination--other bonds are fused at both endpoints. All unselected
 * source bonds stay at source; no line order, graph distance or adjacency between
 * source and destination is required. No singular values are discarded.
 *
 * Collective on comm: graph, selected endpoints and owners must agree on ranks.
 * source, destination and every other must be distinct; the selected bonds must
 * exist and the graph must be simple. bond_neighbors must be nonempty/unique.
 * Tensors and owner entries are required for source, destination and only those
 * other sites with an existing destination--other edge (their axes are fused).
 * Unaffected graph sites need not have tensors/owners. Supplied local tensors must
 * belong to their declared owner and axes must follow incident-edge order; open
 * legs must be represented by graph edges, not implicit extra tensor axes.
 *
 * Labels/owners/local site order stay unchanged. Tensor values and graph/axes
 * are updated together. Undirected edge endpoints may be normalized. Structural
 * errors are rejected before mutation. Backend errors during final absorption
 * may leave partial output, as for graph_tensor_contraction.
 */
template<class TenT, class IntT>
void graph_shift_bonds_to_neighbor(context_handle_t<TenT>& ctx,
    std::vector<std::pair<IntT,IntT>>& edges,std::vector<TenT>& tensors,
    std::vector<IntT>& site_idx,std::map<IntT,IntT>& owners,
    IntT source,IntT destination,const std::vector<IntT>& bond_neighbors,
    MPI_Comm comm) {
    using Index=bond_idx_t<TenT>;
    using Label=bond_label_t<TenT>;
    using RealTensor=real_ten_t<TenT>;
    std::set<IntT> selected;
    std::vector<IntT> required{source,destination};
    if(source==destination || bond_neighbors.empty())
      throw std::invalid_argument("graph bond shift: distinct endpoints and at least one bond required");
    for(auto other:bond_neighbors) {
      if(other==source || other==destination || !selected.insert(other).second ||
         find_edge(edges,source,other)==edges.end())
        throw std::invalid_argument("graph bond shift: invalid selected bond");
      if(find_edge(edges,destination,other)!=edges.end()) required.push_back(other);
    }
    const auto local=detail::validate_graph_tensors(ctx,required,edges,tensors,site_idx,owners,comm);
    int rank;MPI_Comm_rank(comm,&rank);
    int any_invalid=0;
    std::set<IntT> occupied;
    for(auto [a,b]:edges) {occupied.insert(a);occupied.insert(b);}
    for(const auto& entry:owners) occupied.insert(entry.first);
    // Reserve a collision-free temporary label, including labels absent from edges.
    IntT temporary=std::numeric_limits<IntT>::max();
    while(occupied.count(temporary)) {
      if(temporary==std::numeric_limits<IntT>::lowest())
        throw std::invalid_argument("graph bond shift: no temporary site label available");
      --temporary;
    }
    const auto right=source,left=destination;
    const auto bonds=GetSurroundingBondIndex(right,edges);
    std::vector<Index> moving,kept;
    for(std::size_t k=0;k<bonds.size();++k) {
      const auto [a,b]=edges[bonds[k]];
      const auto other=(a==right) ? b : a;
      if(other==left || selected.count(other)) moving.push_back(static_cast<Index>(k));
      else kept.push_back(static_cast<Index>(k));
    }
    const bool outer_product=find_edge(edges,right,left)==edges.end();
    TenT factor,retained,left_expanded;
    int failed=0;
    try {
      if(rank==owners.at(right)) {
        auto t=tcapi::copy(ctx,tensors[local.at(right)]);
        tcapi::List<Index> order(moving.begin(),moving.end());
        order.insert(order.end(),kept.begin(),kept.end());
        tcapi::transpose(ctx,t,order);
        // SVD requires a nonempty column partition, even for a scalar remainder.
        if(kept.empty()) {
          auto shape=tcapi::shape(ctx,t);shape.push_back(1);tcapi::reshape(ctx,t,shape);
        }
        TenT u,vd,s;
        RealTensor sigma;
        tcapi::svd(ctx,t,static_cast<order_t<TenT>>(moving.size()),u,sigma,vd);
        if constexpr(std::is_same_v<TenT,RealTensor>) {
          s=tcapi::copy(ctx,sigma);
        } else {
          s=tcapi::to_cplx(ctx,sigma);
        }
        tcapi::List<Label> lu(moving.size()+1),lf(moving.size()+1);
        std::iota(lu.begin(),lu.end(),0);std::iota(lf.begin(),lf.end(),0);
        lu.back()=-1;
        tcapi::List<Label> ls{-1,static_cast<Label>(moving.size())};
        tcapi::contract(ctx,u,lu,s,ls,factor,lf);
        if(kept.empty()) {
          auto shape=tcapi::shape(ctx,vd);shape.pop_back();tcapi::reshape(ctx,vd,shape);
        }
        tcapi::List<Index> restore(kept.size()+1);
        for(std::size_t k=0;k<kept.size();++k) restore[k]=static_cast<Index>(k+1);
        restore.back()=0;tcapi::transpose(ctx,vd,restore);
        retained=tcapi::move(ctx,vd);
        if(outer_product) {
          auto shape=tcapi::shape(ctx,factor);shape.push_back(1);tcapi::reshape(ctx,factor,shape);
        }
      }
      if(outer_product && rank==owners.at(left)) {
        left_expanded=tcapi::copy(ctx,tensors[local.at(left)]);
        auto shape=tcapi::shape(ctx,left_expanded);shape.push_back(1);
        tcapi::reshape(ctx,left_expanded,shape);
      }
    } catch(const std::exception&) {failed=1;}
    MPI_Allreduce(&failed,&any_invalid,1,MPI_INT,MPI_MAX,comm);
    if(any_invalid) throw std::runtime_error("graph bond shift: SVD/factor preparation failed");
    if(rank==owners.at(right)) {
      tensors[local.at(right)]=tcapi::move(ctx,retained);
      tensors.push_back(tcapi::move(ctx,factor));site_idx.push_back(temporary);
    }
    if(outer_product && rank==owners.at(left))
      tensors[local.at(left)]=tcapi::move(ctx,left_expanded);
    edges=detail::make_graph_shift_edges(edges,right,left,selected,temporary);
    owners.emplace(temporary,owners.at(right));
    graph_tensor_contraction(ctx,edges,tensors,site_idx,owners,left,temporary,left,comm);
}

/**
 * @brief Move one source--other bond to destination--other.
 * @see graph_shift_bonds_to_neighbor for tensor/ownership requirements.
 * The new source--destination bond is always present, even if it did not exist.
 */
template<class TenT, class IntT>
void graph_shift_bond_to_neighbor(context_handle_t<TenT>& ctx,
    std::vector<std::pair<IntT,IntT>>& edges,std::vector<TenT>& tensors,
    std::vector<IntT>& site_idx,std::map<IntT,IntT>& owners,
    IntT source,IntT destination,IntT other,MPI_Comm comm) {
    graph_shift_bonds_to_neighbor(ctx,edges,tensors,site_idx,owners,
                                 source,destination,std::vector<IntT>{other},comm);
}

/** Graph-only counterpart of graph_shift_bonds_to_neighbor.
 * Preserves the exact edge-vector order of the tensor operation, including
 * replacement of the old source--destination edge by an appended SVD edge.
 */
template<class IntT>
void graph_shift_bonds_to_neighbor_for_edges(
    std::vector<std::pair<IntT,IntT>>& edges, IntT source, IntT destination,
    const std::vector<IntT>& others) {
    std::set<IntT> selected;
    if(source==destination || others.empty())
        throw std::invalid_argument("graph bond shift: invalid endpoints or empty selection");
    for(auto other:others)
        if(other==source || other==destination || !selected.insert(other).second ||
           find_edge(edges,source,other)==edges.end())
            throw std::invalid_argument("graph bond shift: invalid selected bond");
    // Same temporary-site construction and absorption as the tensor path.
    std::set<IntT> occupied{source,destination};
    for(auto [a,b]:edges) {occupied.insert(a);occupied.insert(b);}
    IntT temporary=std::numeric_limits<IntT>::max();
    while(occupied.count(temporary)) {
        if(temporary==std::numeric_limits<IntT>::lowest())
            throw std::invalid_argument("graph bond shift: no temporary site label available");
        --temporary;
    }
    auto updated=detail::make_graph_shift_edges(edges,source,destination,selected,temporary);
    graph_tensor_contraction_for_edges(updated,destination,temporary,destination);
    edges=std::move(updated);
}

namespace detail {
  template <typename TenT, typename IntT, typename LevelT>
  void graph_tensor_network_contraction_impl(context_handle_t<TenT> & ctx,
					std::vector<std::pair<IntT,IntT>> & edges,
					std::vector<TenT> & W,
					std::vector<IntT> & site_idx,
					std::map<IntT,IntT> & site_to_mpi_rank,
					const std::vector<std::pair<IntT,IntT>> & edges_schedule,
					const std::map<IntT,LevelT> & site_priority,
					const std::map<IntT,IntT> & label_map,
					MPI_Comm comm, GraphContractionMode mode) {

#ifdef GCHELPER_DEBUG
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);

    if( mpi_rank == 0 ) {
      std::cout << " graph_tensor_network_contraction: input edges =";
      for(const auto & [i,j] : edges) {
	std::cout << " (" << i << "," << j << ")";
      }
      std::cout << std::endl;
      std::cout << " graph_tensor_network_contraction: order of edges to be contracted =";
      std::size_t step = 0;
      for(const auto & [i,j] : edges_schedule) {
	std::cout << ((step == 0) ? " (" : " -> (") << i << "," << j << ")";
	step++;
      }
      std::cout << std::endl;
    }

#endif

    std::size_t num_cont_edges = edges_schedule.size();
    auto edges_schedule_copy = edges_schedule;
    auto site_priority_copy=graph_contraction_priorities(edges_schedule_copy,site_priority,label_map);

#ifdef GCHELPER_DEBUG
    std::size_t step=0;
#endif
    while( !edges_schedule_copy.empty() ) {
      const auto sites=graph_contraction_sites(edges_schedule_copy.front(),site_priority_copy,label_map);
      const auto site_a=sites[0],site_b=sites[1],site_c=sites[2];
#ifdef GCHELPER_DEBUG
      if( mpi_rank == 0 ) {
	std::cout << " graph_tensor_network_contraction: start contraction, site_a, site_b, site_c = "
		  << site_a << ", " << site_b << ", " << site_c << std::endl;
      }
#endif
      graph_tensor_contraction_impl(ctx,edges,W,site_idx,site_to_mpi_rank,
                                    site_a,site_b,site_c,comm,mode);
      graph_tensor_contraction_for_edges(edges_schedule_copy,
					 site_a,site_b,site_c);
#ifdef GCHELPER_DEBUG
      if( mpi_rank == 0 ) {
	std::cout << " graph_tensor_network_contraction: current edges =";
	for(const auto & [i,j] : edges) {
	  std::cout << " (" << i << "," << j << ")";
	}
	std::cout << " at step " << step << std::endl;
	std::cout << " graph_tensor_network_contraction: scheduled edges =";
	if( edges_schedule_copy.empty() ) {
	  std::cout << " none at step " << step << std::endl;
	} else {
	  std::size_t remaining_step=0;
	  for(const auto & [i,j] : edges_schedule_copy) {
	    std::cout << ( (remaining_step==0) ? " (" : " -> (")
		      << i << "," << j << ")";
	    remaining_step++;
	  }
	  std::cout << " at step " << step << std::endl;
	}
	step++;
      }
#endif

    }
  }


} // namespace detail

/**
 * @brief Performs scheduled contractions of a distributed tensor network
 *        represented by a graph.
 *
 * This function contracts the edges specified by @p edges_schedule in order.
 * The tensor network is represented by the graph @p edges, where each vertex
 * corresponds to a site label and each edge represents a contraction between
 * two site tensors.
 *
 * The site tensor with label @c site_label is stored as
 * @c W[adrs] on the MPI rank specified by
 * @c site_to_mpi_rank.at(site_label), where
 * @c site_label = site_idx[adrs]. The graph structure @p edges and the
 * ownership map @p site_to_mpi_rank are assumed to be shared consistently
 * among all MPI ranks.
 *
 * For each edge in @p edges_schedule, the function determines which MPI rank
 * should own the contracted tensor according to @p site_priority. If the
 * priority of one site is higher than that of the other site, the tensor data
 * is transferred to the rank owning the higher-priority site, and the
 * contraction is performed there. Sites not listed in @p site_priority are
 * treated as having the lowest priority.
 *
 * After a contraction between two sites, the resulting site label is updated
 * according to @p label_map. If a contracted site label appears in
 * @p label_map, the corresponding mapped label is used as the new site label.
 * The site labels in @p site_idx, the graph @p edges, and the ownership map
 * @p site_to_mpi_rank are updated consistently to reflect the new graph
 * structure. If no mapping is specified, the original label is kept.
 *
 * @tparam TenT    Tensor type.
 * @tparam IntT    Integer type used for site labels, MPI ranks, and indices.
 * @tparam LevelT  Priority level type.
 *
 * @param ctx
 *     Context handle for the tensor computation interface.
 *
 * @param edges
 *     Graph edges defining the tensor network. This graph is updated as
 *     contractions proceed. The same graph data is assumed to be shared by
 *     all MPI ranks.
 *
 * @param W
 *     Local tensor storage. Each entry @c W[adrs] corresponds to the tensor
 *     associated with the site label @c site_idx[adrs] owned by the local MPI
 *     rank.
 *
 * @param site_idx
 *     Site labels of the tensors stored locally in @p W.
 *
 * @param site_to_mpi_rank
 *     Map from site labels to the MPI ranks that own the corresponding
 *     tensors. This map is updated as contractions proceed and is assumed to
 *     be shared consistently by all MPI ranks.
 *
 * @param edges_schedule
 *     Ordered list of edges to be contracted. This is typically a subset of
 *     @p edges. The contractions are performed in the order given by this
 *     vector.
 *
 * @param site_priority
 *     Priority map used to determine the destination rank of each contraction.
 *     When contracting sites @c site_a and @c site_b, if
 *     @c site_priority.at(site_a) > site_priority.at(site_b), the tensor on
 *     @c site_b is sent to the rank owning @c site_a and the contraction is
 *     performed there. Sites not present in this map are treated as having the
 *     lowest priority.
 *
 * @param label_map
 *     Map specifying the site label of the tensor resulting from a contraction.
 *     If a site label is present in this map, the contracted site is relabeled
 *     to the mapped value. If no entry is found, the original label is kept.
 *
 * @param comm
 *     MPI communicator over which tensor transfers and contractions are
 *     coordinated.
 */
template<class TenT, class IntT, class LevelT>
void graph_tensor_network_contraction(context_handle_t<TenT>& ctx,
    std::vector<std::pair<IntT,IntT>>& edges, std::vector<TenT>& tensors,
    std::vector<IntT>& site_idx, std::map<IntT,IntT>& owners,
    const std::vector<std::pair<IntT,IntT>>& schedule,
    const std::map<IntT,LevelT>& priority, const std::map<IntT,IntT>& labels,
    MPI_Comm comm) {
    detail::graph_tensor_network_contraction_impl(ctx,edges,tensors,site_idx,owners,
        schedule,priority,labels,comm,detail::GraphContractionMode::complete);
}

namespace detail {
// Sampling-only path: absent operands/boundaries are explicitly owner=-1.
// Existing common-neighbor tensors still undergo exactly the same bond fusion.
template<class TenT, class IntT, class LevelT>
void graph_tensor_network_partial_contraction(context_handle_t<TenT>& ctx,
    std::vector<std::pair<IntT,IntT>>& edges, std::vector<TenT>& tensors,
    std::vector<IntT>& site_idx, std::map<IntT,IntT>& owners,
    const std::vector<std::pair<IntT,IntT>>& schedule,
    const std::map<IntT,LevelT>& priority, const std::map<IntT,IntT>& labels,
    MPI_Comm comm) {
    graph_tensor_network_contraction_impl(ctx,edges,tensors,site_idx,owners,
        schedule,priority,labels,comm,GraphContractionMode::partial_environment);
}
} // namespace detail

}

#endif
