/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/ptns/function/bmps.h
@brief Boundary-MPS grouping, absorption, bond relocation and BP truncation.
@details Graph edges are ordered: tensor axes follow the order of incident edges.
Planning and graph-only operations share their selection rules with tensor
operations. Absorption and bond relocation do not truncate; call
bmps_bp_truncation() separately to limit boundary bond dimensions.
*/
#ifndef TNBP_PTNS_FUNCTION_BMPS_H
#define TNBP_PTNS_FUNCTION_BMPS_H

#include "tnbp/ptns/function/bp.h"
#include "tnbp/ptns/function/init.h"
#include "tnbp/ptns/function/truncation.h"
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <set>
#include <queue>
#include <tuple>
#include <functional>
#include <stdexcept>
#include <ostream>
#include <utility>
#include <vector>

#include "tnbp/framework/typedef.h"
#include "tnbp/framework/graph.h"
#include "tnbp/framework/graph_contraction.h"
#include "tnbp/framework/mpiutility.h"

namespace tnbp {

  /**
   * @brief Build successive graph layers from an initial boundary.
   * @tparam IntT Site-label integer type.
   * @param[in] edges Undirected graph edges.
   * @param[in] initial_boundary First layer, in the requested boundary order.
   * @param[out] lines Replaced with the layers returned by BuildContractionLayers().
   *
   * The first layer preserves initial_boundary; subsequent layers are sorted by
   * site label. Components unreachable from the initial boundary are not included.
   * This wrapper adds no input validation to BuildContractionLayers().
   */
  template <typename IntT>
  void bmps_setup_lines(const std::vector<std::pair<IntT,IntT>> & edges,
			const std::vector<IntT> & initial_boundary,
			std::vector<std::vector<IntT>> & lines) {
    lines = BuildContractionLayers(edges,initial_boundary);
  }

  /**
   * @brief Partition one neighboring line among ordered target sites.
   * @tparam IntT Site-label integer type.
   * @param[in] edges Original graph edges used to select absorption paths.
   * @param[in] targets Ordered anchors of the destination line.
   * @param[in] neighbors Sites of one neighboring line, in tie-breaking order.
   * @return One group per target, starting with that target, followed by its
   * assigned neighbors in absorption order. Each neighbor occurs exactly once.
   * @throws std::invalid_argument If the lists overlap, contain duplicates, or
   * some neighbor cannot be reached from a target through the neighboring line.
   *
   * In target order, reserve the first still-unassigned direct neighbor for each
   * anchor, if available. Then grow a deterministic breadth-first forest, ordering
   * candidates by path length, target position, and neighbor position. Reservations
   * are retained; this is not an unconstrained nearest-anchor partition.
   * Paths may use their own anchor and neighboring-line sites only, never another
   * target as an intermediate. The input graph and lists are not modified.
   * @see bmps_absorption_schedule
   */
  template<class IntT>
  std::vector<std::vector<IntT>> bmps_neighbor_groups(
      const std::vector<std::pair<IntT,IntT>>& edges,
      const std::vector<IntT>& targets,const std::vector<IntT>& neighbors) {
    std::set<IntT> target_set(targets.begin(),targets.end());
    std::map<IntT,std::size_t> index;
    for(std::size_t k=0;k<neighbors.size();++k)
      if(target_set.count(neighbors[k]) || !index.emplace(neighbors[k],k).second)
        throw std::invalid_argument("automatic grouping: overlapping or duplicate sites");
    if(target_set.size()!=targets.size())
      throw std::invalid_argument("automatic grouping: duplicate target");
    std::vector<std::vector<IntT>> groups;
    for(auto site:targets) groups.push_back({site});
    std::vector<bool> assigned(neighbors.size(),false);
    using Candidate=std::tuple<std::size_t,std::size_t,std::size_t>; // distance, anchor, site position
    std::priority_queue<Candidate,std::vector<Candidate>,std::greater<Candidate>> frontier;
    auto offer=[&](IntT site,std::size_t anchor,std::size_t distance) {
      for(auto [a,b]:edges) {
        if(a!=site && b!=site) continue;
        const auto it=index.find(a==site ? b : a);
        if(it!=index.end() && !assigned[it->second])
          frontier.emplace(distance,anchor,it->second);
      }
    };
    for(std::size_t j=0;j<targets.size();++j) {
      for(std::size_t k=0;k<neighbors.size();++k)
        if(!assigned[k] && find_edge(edges,targets[j],neighbors[k])!=edges.end()) {
          assigned[k]=true;groups[j].push_back(neighbors[k]);
          break;
        }
    }
    for(std::size_t j=0;j<targets.size();++j) {
      offer(targets[j],j,1);
      if(groups[j].size()>1) offer(groups[j][1],j,2);
    }
    while(!frontier.empty()) {
      const auto [distance,j,k]=frontier.top();frontier.pop();
      if(assigned[k]) continue;
      assigned[k]=true;groups[j].push_back(neighbors[k]);
      offer(neighbors[k],j,distance+1);
    }
    if(std::find(assigned.begin(),assigned.end(),false)!=assigned.end())
      throw std::invalid_argument("automatic grouping: neighboring component has no path to a target");
    return groups;
  }

  /**
   * @brief Validate and plan absorption of an input line into an output line.
   * @tparam IntT Site-label integer type.
   * @param[in] input Sites to absorb, each exactly once.
   * @param[in] output Ordered destination sites, disjoint from input.
   * @param[in] edges Current ordered graph edges; never modified by this function.
   * @param[in] groups Optional groups, one per output site and starting with that
   * site. A null pointer selects bmps_neighbor_groups(edges, output, input).
   * @return Ordered pairs (destination site, absorbed site).
   * @throws std::invalid_argument If lines or groups are invalid, an input site
   * is omitted, or a scheduled absorption has no edge in the evolving graph.
   *
   * With explicit groups, iterate groups and their members in supplied order.
   * Members outside input, such as the opposite neighboring line in a saved
   * three-line partition, are validated but not absorbed. Members must be unique
   * across groups and present in edges; other output sites cannot be members.
   * All input sites must be covered. A nonempty input requires a nonempty output.
   *
   * Each planned contraction is simulated on a graph copy. Thus planning errors
   * leave the caller's graph untouched; tensor shapes and owners are not checked.
   * Sampling supplies its saved three-line groups, avoiding directional regrouping.
   * @see bmps_increment_line
   * @see bmps_increment_line_for_edges
   */
  template<class IntT>
  std::vector<std::pair<IntT,IntT>> bmps_absorption_schedule(
      const std::vector<IntT>& input, const std::vector<IntT>& output,
      const std::vector<std::pair<IntT,IntT>>& edges,
      const std::vector<std::vector<IntT>>* groups=nullptr) {
    const std::set<IntT> inputs(input.begin(),input.end());
    const std::set<IntT> outputs(output.begin(),output.end());
    if(inputs.size()!=input.size() || outputs.size()!=output.size())
      throw std::invalid_argument("BMPS: duplicate line site");
    for(auto site:input) if(outputs.count(site))
      throw std::invalid_argument("BMPS: input and output lines overlap");
    if(!input.empty() && output.empty())
      throw std::invalid_argument("BMPS: output line is empty");
    auto temporary=edges;
    auto remaining=inputs;
    std::vector<std::pair<IntT,IntT>> schedule;
    auto absorb=[&](IntT target,IntT site) {
      if(!remaining.count(site) || find_edge(temporary,target,site)==temporary.end())
        throw std::invalid_argument("BMPS: absorption requires an available site and edge");
      schedule.emplace_back(target,site);
      graph_tensor_contraction_for_edges(temporary,target,site,target);
      remaining.erase(site);
    };
    if(groups) {
      if(groups->size()!=output.size())
        throw std::invalid_argument("BMPS: group count differs from target line");
      std::set<IntT> members;
      const auto sites=GetSiteIndexFromBond(edges);
      for(std::size_t j=0;j<output.size();++j) {
        const auto& group=(*groups)[j];
        if(group.empty() || group.front()!=output[j])
          throw std::invalid_argument("BMPS: group must start with its target site");
        for(std::size_t k=0;k<group.size();++k) {
          const auto site=group[k];
          if(!members.insert(site).second ||
             std::find(sites.begin(),sites.end(),site)==sites.end() ||
             (k && outputs.count(site)))
            throw std::invalid_argument("BMPS: duplicate, unknown or misplaced group site");
          // Groups may also contain the opposite neighboring line. Those sites
          // remain untouched during this direction of boundary absorption.
          if(k && inputs.count(site)) absorb(output[j],site);
        }
      }
    } else {
      const auto automatic=bmps_neighbor_groups(edges,output,input);
      for(std::size_t j=0;j<output.size();++j)
        for(std::size_t k=1;k<automatic[j].size();++k)
          absorb(output[j],automatic[j][k]);
    }
    if(!remaining.empty()) throw std::invalid_argument("BMPS: grouping omits an input-line site");
    return schedule;
  }

  /**
   * @brief Select long boundary bonds in a shared right-to-left sweep.
   * @tparam IntT Site-label integer type.
   * @tparam Shift Callable accepting (source, destination, bond_neighbors).
   * @param[in] line Boundary sites in left-to-right order.
   * @param[in] edges Graph observed at each sweep step. The callback must update
   * this same graph through a mutable reference if subsequent steps depend on it.
   * @param[in] shift Operation moving the selected source bonds to its left neighbor.
   * @pre line contains unique site labels and edges describes a valid simple graph.
   *
   * At each source, select endpoints strictly left of its immediate left neighbor,
   * in incident-edge order. Invoke shift only when this list is nonempty. Edges
   * outside the line are never selected. Empty and singleton lines do no work.
   * This helper performs selection only; it does not validate inputs or access
   * MPI/tensors. Tensor and graph-only wrappers supply the actual mutation.
   * @see bmps_shift_bonds_to_neighbors
   * @see bmps_shift_bonds_to_neighbors_for_edges
   */
  template<class IntT, class Shift>
  void bmps_shift_bond_sweep(const std::vector<IntT>& line,
      const std::vector<std::pair<IntT,IntT>>& edges, Shift shift) {
    std::map<IntT,std::size_t> position;
    for(std::size_t j=0;j<line.size();++j) position.emplace(line[j],j);
    for(std::size_t j=line.size();j>1;--j) {
      const auto right=line[j-1],left=line[j-2];
      std::vector<IntT> targets;
      for(auto e:GetSurroundingBondIndex(right,edges)) {
        const auto [a,b]=edges[e];
        const auto other=(a==right) ? b : a;
        const auto p=position.find(other);
        if(p!=position.end() && p->second<j-2) targets.push_back(other);
      }
      if(!targets.empty())
        shift(right,left,targets);
    }
  }

  /**
   * @brief Move long boundary bonds onto adjacent sites without truncation.
   * @tparam TenT TCAPI tensor type.
   * @tparam IntT Site-label and owner-rank integer type.
   * @param[in,out] ctx Tensor backend context.
   * @param[in] line Boundary sites in left-to-right order.
   * @param[in,out] edges Ordered graph edges, updated together with tensor axes.
   * @param[in,out] tensors Local tensors; only boundary tensors are changed.
   * @param[in,out] site_idx Local tensor labels; values and order are preserved.
   * @param[in,out] owners Site-to-MPI-rank map; ownership is preserved.
   * @param[in] comm MPI communicator; all ranks participate.
   * @pre edges, line and owners agree across ranks. The graph is simple and each
   * line tensor is supplied on its declared owner. Axes follow incident-edge order;
   * external legs are graph edges, not implicit extra axes.
   * @throws std::invalid_argument If structural or tensor-shape validation fails.
   *
   * Sweep right to left using bmps_shift_bond_sweep(). For each site with a bond
   * skipping its left neighbor, use graph_shift_bonds_to_neighbor() to move the
   * selected bonds and any existing source-destination bond into the left factor.
   * A full TCAPI SVD retains V-dagger at the source and absorbs U*S at the left
   * neighbor. An absent neighbor bond uses a temporary dimension-one bond; common
   * neighbor bonds are fused. External bonds stay at their original sites.
   *
   * Only line tensors and their owner entries are required, including common
   * neighbors whose axes may be fused. Unaffected tensors may be omitted.
   * Already adjacent bonds alone do not trigger SVD. No singular values are
   * discarded: this is not a full MPS canonicalization, and dimensions may grow.
   * Structural/shape validation precedes the sweep collectively. Later backend
   * failures may leave partial results; the operation is not transactional.
   * @see bmps_shift_bonds_to_neighbors_for_edges
   * @see bmps_bp_truncation
   */
  template<class TenT, class IntT>
  void bmps_shift_bonds_to_neighbors(context_handle_t<TenT>& ctx,
      const std::vector<IntT>& line, std::vector<std::pair<IntT,IntT>>& edges,
      std::vector<TenT>& tensors, std::vector<IntT>& site_idx,
      std::map<IntT,IntT>& owners, MPI_Comm comm) {
    detail::validate_graph_tensors(ctx,line,edges,tensors,site_idx,owners,comm);
    bmps_shift_bond_sweep(line,edges,[&](IntT source,IntT destination,const std::vector<IntT>& targets) {
      graph_shift_bonds_to_neighbor(ctx,edges,tensors,site_idx,owners,source,destination,targets,comm);
    });
  }

  /**
   * @brief Apply boundary bond relocation to the graph alone.
   * @tparam IntT Site-label integer type.
   * @param[in] line Boundary sites in left-to-right order; labels must be unique.
   * @param[in,out] edges Ordered simple graph edges.
   *
   * Uses the same evolving-graph selection and edge updates as the tensor version,
   * without tensors, owners, SVD or MPI. Valid inputs produce the same ordered edge
   * list. This wrapper does not perform the tensor version's initial validation;
   * individual graph shifts validate their own inputs. A later failure can leave
   * earlier shifts applied. Empty and singleton lines perform no shifts.
   * @see bmps_shift_bonds_to_neighbors
   */
  template<class IntT>
  void bmps_shift_bonds_to_neighbors_for_edges(const std::vector<IntT>& line,
      std::vector<std::pair<IntT,IntT>>& edges) {
    bmps_shift_bond_sweep(line,edges,[&](IntT source,IntT destination,const std::vector<IntT>& targets) {
      graph_shift_bonds_to_neighbor_for_edges(edges,source,destination,targets);
    });
  }


  /**
   * @brief Absorb a boundary into the next line and remove long line bonds.
   * @tparam TenT TCAPI tensor type.
   * @tparam IntT Site-label and owner-rank integer type.
   * @param[in,out] ctx Tensor backend context.
   * @param[in] line_in Boundary sites to absorb and remove.
   * @param[in] line_out Ordered surviving destination sites.
   * @param[in,out] edges Ordered graph edges, including relevant external edges.
   * @param[in,out] tensors Local tensors updated by absorption and bond relocation.
   * @param[in,out] site_idx Local labels updated when input sites are removed.
   * @param[in,out] owners Site-to-MPI-rank map updated when sites are absorbed.
   * Each surviving destination keeps its original owner.
   * @param[in] comm MPI communicator; all ranks participate.
   * @param[in] groups Optional saved groups; see bmps_absorption_schedule().
   * @pre Graph, lines, groups and owner metadata agree across ranks. Supply both
   * operands and every common neighbor affected by axis fusion at each contraction,
   * plus all output-line tensors, on their declared owners. Unaffected tensors may
   * be omitted even when their edges appear in edges. Axes follow incident-edge order.
   *
   * First plan on a graph copy, then contract each (target, site) into target using
   * graph_tensor_contraction(). Finally call bmps_shift_bonds_to_neighbors() on
   * line_out. Sites from the opposite neighboring line in groups are not absorbed.
   * No BP or dimension truncation is performed here.
   * @throws std::invalid_argument If planning or a delegated operation rejects input.
   * @note Planning errors precede all mutation; tensor validation occurs in the
   * individual operations. A failure after absorption starts may leave partial output.
   * @see bmps_increment_line_for_edges
   * @see bmps_bp_truncation
   */
  template <typename TenT, typename IntT>
  void bmps_increment_line(context_handle_t<TenT>& ctx,
      const std::vector<IntT>& line_in, const std::vector<IntT>& line_out,
      std::vector<std::pair<IntT,IntT>>& edges, std::vector<TenT>& tensors,
      std::vector<IntT>& site_idx, std::map<IntT,IntT>& owners, MPI_Comm comm,
      const std::vector<std::vector<IntT>>* groups=nullptr) {
    const auto schedule=bmps_absorption_schedule(line_in,line_out,edges,groups);
#ifdef BMPS_DEBUG
    int rank;MPI_Comm_rank(comm,&rank);
#endif
    for(const auto& [target,site]:schedule) {
#ifdef BMPS_DEBUG
      if(rank==0) std::cout << "BMPS absorption: " << site << " -> " << target << std::endl;
#endif
      graph_tensor_contraction(ctx,edges,tensors,site_idx,owners,target,site,target,comm);
    }
    bmps_shift_bonds_to_neighbors(ctx,line_out,edges,tensors,site_idx,owners,comm);
  }

  /**
   * @brief Simulate boundary absorption and bond relocation on ordered edges.
   * @tparam IntT Site-label integer type.
   * @param[in] line_in Sites to absorb and remove.
   * @param[in] line_out Ordered surviving destination sites.
   * @param[in,out] edges Ordered simple graph edges.
   * @param[in] groups Optional saved groups; see bmps_absorption_schedule().
   *
   * Shares the absorption schedule and right-to-left bond selection with
   * bmps_increment_line(), producing the same edge order for valid inputs.
   * No tensors, owner metadata, MPI operations or truncation are involved.
   * @throws std::invalid_argument If planning or a delegated graph operation fails.
   * @note Planning uses a copy, but execution is not transactional: a later error
   * may leave earlier edge updates applied.
   * @see bmps_increment_line
   */
  template <typename IntT>
  void bmps_increment_line_for_edges(const std::vector<IntT>& line_in,
      const std::vector<IntT>& line_out, std::vector<std::pair<IntT,IntT>>& edges,
      const std::vector<std::vector<IntT>>* groups=nullptr) {
    const auto schedule=bmps_absorption_schedule(line_in,line_out,edges,groups);
    for(const auto& [target,site]:schedule)
      graph_tensor_contraction_for_edges(edges,target,site,target);
    bmps_shift_bonds_to_neighbors_for_edges(line_out,edges);
  }

  /**
   * @brief Run belief propagation and truncate bonds within one boundary line.
   * @tparam TenT TCAPI tensor type.
   * @tparam IntT Site-label and owner-rank integer type.
   * @tparam StepT Integer type for the BP iteration count.
   * @param[in] ctx Tensor backend context handle.
   * @param[in] edges Full ordered graph edges; connectivity and order are preserved.
   * @param[in] lines Sites of one boundary line (not a collection of lines).
   * @param[in,out] T Local tensors, updated by truncation.
   * @param[in] site_idx Labels corresponding to local entries of T.
   * @param[in] site_to_mpi_rank Global site-to-MPI-rank map; ownership is unchanged.
   * @param[in,out] E BP messenger tensors; initialized when do_init is true,
   * otherwise reused as the initial messages. Updated by BP and truncation.
   * @param[in,out] edge_idx Messenger edge indices relative to the induced line
   * edge list; initialized with E when do_init is true.
   * @param[in] max_bp_steps Maximum BP iterations before truncation.
   * @param[in] bp_tolerance Stop BP when the mean normalized, trace-gauged successive-message
   * distance is strictly below this value.
   * @param[in] max_bond_dim Maximum retained bond dimension passed to Truncation().
   * @param[in] sv_min Singular-value threshold passed to Truncation().
   * @param[in] tg_err Target truncation error passed to Truncation().
   * @param[out] res_bond_dim Retained dimensions reported by Truncation().
   * @param[out] res_trunc_err Truncation errors reported by Truncation().
   * @param[in] comm MPI communicator; all ranks participate.
   * @param[in] do_init Whether to initialize messenger tensors (default true).
   * @pre Graph, line and owner metadata agree across ranks. Each line tensor is
   * available on its declared owner and axes follow the full incident-edge order.
   * When reusing E/edge_idx, they must match the induced graph and current dimensions.
   *
   * Extract the edges with both endpoints in lines. Temporarily transpose each
   * line tensor so these internal axes precede external axes, run BP with local
   * message updates followed by rank-pair boundary exchanges, then call
   * Truncation() on the induced graph. On normal return,
   * restore tensor axes to the original incident-edge order. External axes are
   * not truncated. Reaching max_bp_steps still proceeds to truncation; convergence
   * is not required, and zero steps skips BP updates.
   * @note This routine assumes valid inputs rather than performing the shared graph
   * validation. An exception after transposition may leave tensors in that temporary
   * axis order. Site labels and owners are never changed.
   * @see bmps_increment_line
   */
  template <typename TenT, typename IntT, typename StepT>
  void bmps_bp_truncation(context_handle_t<TenT>& ctx,
			  const std::vector<std::pair<IntT,IntT>> & edges,
			  const std::vector<IntT> & lines,
			  std::vector<TenT> & T,
			  const std::vector<IntT> & site_idx,
			  const std::map<IntT,IntT> & site_to_mpi_rank,
			  std::vector<TenT> & E,
			  std::vector<IntT> & edge_idx,
			  StepT max_bp_steps,
			  real_t<TenT> bp_tolerance,
			  bond_dim_t<TenT> max_bond_dim,
			  real_t<TenT> sv_min,
			  real_t<TenT> tg_err,
			  std::vector<bond_dim_t<TenT>> & res_bond_dim,
			  std::vector<real_t<TenT>> & res_trunc_err,
			  MPI_Comm comm,
			  bool do_init = true) {
    
    using RealT = typename tcapi::tensor_traits<TenT>::real_t;
    using BondDimT = typename tcapi::tensor_traits<TenT>::bond_dim_t;
    using BondIdxT = typename tcapi::tensor_traits<TenT>::bond_idx_t;
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    auto edges_for_bp = extract_induced_edges(edges,lines);
    // It is necessary to transpose tensors
    std::vector<std::vector<BondIdxT>> trs_order(lines.size());
    for(std::size_t adrs=0; adrs < lines.size(); adrs++) {
      auto site_i = lines[adrs];
      auto mpi_rank_i = site_to_mpi_rank.at(site_i);
      if( mpi_rank == mpi_rank_i ) {
	auto bond_all = GetSurroundingBondIndex(site_i,edges);
	auto bond_mps = GetSurroundingBondIndex(site_i,edges_for_bp);
	trs_order[adrs].resize(bond_all.size());
	BondIdxT kpos = static_cast<BondIdxT>(bond_mps.size());
	for(std::size_t k=0; k < bond_all.size(); k++) {
	  BondIdxT mpos = bond_all.size();
	  for(std::size_t m=0; m < bond_mps.size(); m++) {
	    if( same_edge(edges_for_bp[bond_mps[m]],
			  edges[bond_all[k]]) ) {
	      mpos = m;
	      break;
	    }
	  }
	  if( mpos == bond_all.size() ) {
	    trs_order[adrs][kpos++] = k;
	  } else {
	    trs_order[adrs][mpos] = k;
	  }
	}
	auto it_site_adrs = std::find(site_idx.begin(),
				      site_idx.end(),
				      site_i);
	auto site_adrs = std::distance(site_idx.begin(),
				       it_site_adrs);
	tcapi::transpose(ctx,T[site_adrs],trs_order[adrs]);
      }
    }
    
    if( do_init ) {
      init_edge_messenger_tensors(ctx,edges_for_bp,T,
				  site_idx,site_to_mpi_rank,
				  E,edge_idx,comm);
    }
    int bp_rank; MPI_Comm_rank(comm,&bp_rank);
    const auto bp_order=MakeBPCommunicationOrder(edges_for_bp,site_to_mpi_rank,bp_rank);
    for(StepT step=0; step < max_bp_steps; step++) {
      std::vector<TenT> E_new;
      UpdateBPMessages(ctx,edges_for_bp,T,site_idx,E,edge_idx,comm,E_new);
      ExchangeBPMessages(ctx,edges_for_bp,site_to_mpi_rank,edge_idx,bp_order,comm,E_new);
      RealT res_bp_err;
      BeliefPropagationDifference(ctx,edges_for_bp,site_to_mpi_rank,
                                 E,E_new,edge_idx,comm,res_bp_err);
      E.swap(E_new);
#ifdef BMPS_DEBUG
      if( mpi_rank == 0 ) {
	std::cout << " bp step " << step << ": error = " << res_bp_err << std::endl;
      }
#endif
      if( res_bp_err < bp_tolerance ) {
	break;
      }
    }
    Truncation(ctx,edges_for_bp,T,site_idx,site_to_mpi_rank,
	       E,edge_idx,comm,max_bond_dim,sv_min,tg_err,
	       res_bond_dim,res_trunc_err);

    for(std::size_t adrs=0; adrs < lines.size(); adrs++) {
      auto site_i = lines[adrs];
      auto mpi_rank_i = site_to_mpi_rank.at(site_i);
      if( mpi_rank == mpi_rank_i ) {
	std::vector<BondIdxT> rev_order(trs_order[adrs].size());
	for(std::size_t k=0; k < trs_order[adrs].size(); k++) {
	  rev_order[trs_order[adrs][k]] = k;
	}
	auto it_site_adrs = std::find(site_idx.begin(),
				      site_idx.end(),
				      site_i);
	auto site_adrs = std::distance(site_idx.begin(),
				       it_site_adrs);
	tcapi::transpose(ctx,T[site_adrs],rev_order);
      }
    }
  }
			  
  
  
}

#endif
