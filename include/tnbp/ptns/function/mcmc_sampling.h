/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/ptns/function/mcmc_sampling.h
@brief Boundary-MPS Markov-chain Monte Carlo initialization and sampling.
@details Original tensors W have virtual axes in incident-edge order followed by
one physical axis. Attached tensors V have the physical state fixed and only
virtual axes. Tensor and state arrays are rank-local, indexed by site_idx;
graphs, ordered lines and grouping metadata are shared across MPI ranks.
Upper/lower E arrays contain BP messages, not the left/right sampling environments.
Sampling reuses the grouping fixed during initialization in both sweep directions.
*/
#ifndef TNBP_PTNS_FUNCTION_MCMC_SAMPLING_H
#define TNBP_PTNS_FUNCTION_MCMC_SAMPLING_H

#include "tnbp/ptns/function/mcmc_plan.h"
#include "tnbp/ptns/function/mcmc_heat_bath.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>
#include "tnbp/framework/typedef.h"
#include "tnbp/framework/graph.h"
#include "tnbp/framework/graph_contraction.h"
#include "tnbp/framework/mpiutility.h"

#ifdef MCMC_DEBUG
#include <unistd.h>
#endif

namespace tnbp {

  /**
   * @brief Validate the requested update method.
   * @param[in] method 0 for Metropolis, 1 for Barker, 2 for full-state heat bath.
   *
   * @throws std::invalid_argument If method is not 0, 1 or 2.
   */
  inline void validate_mcmc_sampling_method(int method) {
    if(method!=0 && method!=1 && method!=2)
      throw std::invalid_argument("sampling_method must be 0 (Metropolis), 1 (Barker), or 2 (full-state heat bath)");
  }

  /**
   * @brief Compute the Metropolis acceptance probability.
   * @tparam RealT Real floating-point type.
   * @param[in] current Finite nonnegative absolute amplitude of the current state.
   * @param[in] proposed Finite nonnegative absolute amplitude of the proposed state.
   *
   * @return min(1, proposed^2/current^2), or zero if both amplitudes vanish.
   * Compare amplitudes first and square only a ratio below one to avoid separately
   * squaring very large or small amplitudes.
   * A positive proposal from a zero-weight state is accepted with probability one.
   * @pre Inputs satisfy the finite/nonnegative contract; no validation is performed.
   */
  template <typename RealT>
  RealT mcmc_pair_metropolis_probability(RealT current, RealT proposed) {
    if(proposed==RealT(0)) return RealT(0);
    if(proposed>=current) return RealT(1);
    const RealT ratio=proposed/current;
    return ratio*ratio;
  }

  /**
   * @brief Compute the pairwise heat-bath acceptance probability.
   * @tparam RealT Real floating-point type.
   * @param[in] current Finite nonnegative absolute amplitude of the current state.
   * @param[in] proposed Finite nonnegative absolute amplitude of the proposed state.
   *
   * @return proposed^2/(current^2 + proposed^2), or zero if both vanish.
   * Rescale amplitudes before squaring. For physical dimension two this is a full
   * heat-bath update; for larger dimensions it is a two-state Barker update.
   * A positive proposal from a zero-weight state is accepted with probability one.
   * @pre Inputs satisfy the finite/nonnegative contract; no validation is performed.
   */
  template <typename RealT>
  RealT mcmc_pair_heat_bath_probability(RealT current, RealT proposed) {
    const RealT scale=std::max(current,proposed);
    if(scale==RealT(0)) return RealT(0);
    const RealT a=current/scale, b=proposed/scale;
    return b*b/(a*a+b*b);
  }


  /**
   * @brief Create a one-hot physical-state tensor.
   * @tparam TenT TCAPI tensor type.
   * @tparam IntT Signed integer type for site labels, states and indices.
   * @tparam DimT Integer type for the physical dimension; may differ from IntT.
   * @param[in,out] ctx TCAPI tensor backend context.
   * @param[in] s Selected physical state, with 0 <= s < d.
   * @param[in] d Positive physical dimension.
   * @param[out] U Rank-one tensor of length d, one at s and zero elsewhere.
   *
   * @pre d and s are valid; bounds are not checked. This is a local operation.
   */
  template <typename TenT, typename IntT, typename DimT>
  void state_tensor(context_handle_t<TenT> & ctx,
		    IntT s,
		    DimT d,
		    TenT & U) {
    using ElemT = typename tcapi::tensor_traits<TenT>::elem_t;
    using BondDimT = typename tcapi::tensor_traits<TenT>::bond_dim_t;
    using ShapeT = typename tcapi::tensor_traits<TenT>::shape_t;
    using CoorsT = typename tcapi::tensor_traits<TenT>::elem_coors_t;
    ShapeT shape_u(1,static_cast<BondDimT>(d));
    List<ElemT> data_u(d,ElemT(0.0));
    data_u[s] = ElemT(1.0);
    auto it_data_u = data_u.begin();
    U = tcapi::assign_from_range<TenT>(ctx,shape_u,it_data_u,
				     [&shape_u](const CoorsT & coors) {
				       return address_from_coor(shape_u,coors);
				     });
  }

  /**
   * @brief Build ordered layers from an initial graph boundary.
   * @tparam IntT Signed integer type for site labels, states and indices.
   * @param[in] edges Ordered undirected graph edges; virtual axes follow incident-edge order.
   * @param[in] initial_boundary First boundary layer, in the requested site order.
   * @param[out] lines Layers returned by bmps_setup_lines().
   *
   * Delegates to bmps_setup_lines(); disconnected unreachable components are not
   * included. This function does not initialize a workspace or choose groups.
   * @see bmps_setup_lines
   */
  template <typename IntT>
  void mcmc_init_line(const std::vector<std::pair<IntT,IntT>> & edges,
		      const std::vector<IntT> & initial_boundary,
		      std::vector<std::vector<IntT>> & lines) {
    bmps_setup_lines(edges,initial_boundary,lines);
  }

  /**
   * @brief Fix each local tensor physical leg to the current sampled state.
   * @tparam TenT TCAPI tensor type.
   * @tparam IntT Signed integer type for site labels, states and indices.
   * @param[in,out] ctx TCAPI tensor backend context.
   * @param[in] edges Graph edges, retained for interface consistency; not read here.
   * @param[in] W Original rank-local tensors, with the physical axis last; never modified.
   * @param[in] site_idx Site labels in the same local order as the tensor/state arrays.
   * @param[in] site_to_mpi_rank Owner map, retained for interface consistency; not read here.
   * @param[in] site_states Physical states indexed by local site_idx position, not by site label.
   * @param[in] lines Ordered nonempty lines, from lower to upper boundary.
   * @param[out] V Attached tensors in W order, with the final physical axis contracted away.
   *
   * Only local sites found in lines are assigned. No MPI communication is performed.
   * @pre W, site_idx and site_states have matching local lengths; lines cover all
   * local sites, each state is in range, and every W tensor has a final physical axis.
   * @see state_tensor
   */
  template <typename TenT, typename IntT>
  void mcmc_init_attachment(context_handle_t<TenT> & ctx,
			    const std::vector<std::pair<IntT,IntT>> & edges,
			    const std::vector<TenT> & W,
			    const std::vector<IntT> & site_idx,
			    const std::map<IntT,IntT> & site_to_mpi_rank,
			    const std::vector<IntT> & site_states,
			    const std::vector<std::vector<IntT>> & lines,
			    std::vector<TenT> & V) {

    using BondLabelT = typename tcapi::tensor_traits<TenT>::bond_label_t;
    V.resize(W.size());
    for(std::size_t l=0; l < lines.size(); l++) {
      for(std::size_t i=0; i < lines[l].size(); i++) {
	auto it_site_address = std::find(site_idx.begin(),site_idx.end(),lines[l][i]);
	if( it_site_address != site_idx.end() ) {
	  auto site_address = std::distance(site_idx.begin(),it_site_address);
	  auto shape_w = tcapi::shape(ctx,W[site_address]);
	  auto order_w = tcapi::order(ctx,W[site_address]);
	  TenT S;
	  auto dim_s = shape_w[order_w-1];
	  state_tensor(ctx,site_states[site_address],dim_s,S);
	  auto shape_s = tcapi::shape(ctx,S);
	  auto order_s = tcapi::order(ctx,S);
	  auto order_v = order_w-1;
	  std::vector<BondLabelT> label_w(order_w);
	  std::vector<BondLabelT> label_s(order_s);
	  std::vector<BondLabelT> label_v(order_v);
	  std::iota(label_w.begin(),label_w.end(),0);
	  label_w[order_w-1] = -1;
	  label_s[0] = -1;
	  std::iota(label_v.begin(),label_v.end(),0);
	  tcapi::contract(ctx,W[site_address],label_w,
			S,label_s,
			V[site_address],label_v);
	}
      }
    }
  }

  /**
   * @brief Replay boundary absorption to obtain a target line graph.
   * @tparam IntT Signed integer type for site labels, states and indices.
   * @param[in] edges Ordered undirected graph edges; virtual axes follow incident-edge order.
   * @param[in] lines Ordered nonempty lines, from lower to upper boundary.
   * @param[in] target_line Index of the target in lines.
   * @param[out] target_edges Graph after outer regions are absorbed into the adjacent boundary lines.
   * @param[in] grouping Optional fixed three-line grouping. Null selects automatic grouping from the original graph.
   *
   * The target and its immediate neighboring lines remain. Complete graph-only BMPS
   * increments include bond relocation, preserving the tensor path edge order.
   * @pre target_line is valid and a supplied grouping matches edges and lines.
   * @see bmps_increment_line_for_edges
   */
  template <typename IntT>
  void mcmc_init_edges(const std::vector<std::pair<IntT,IntT>> & edges,
		       const std::vector<std::vector<IntT>> & lines,
		       IntT target_line,
		       std::vector<std::pair<IntT,IntT>> & target_edges,
                       const MCMCGrouping<IntT>* grouping=nullptr) {

    MCMCGrouping<IntT> automatic;
    if(!grouping) { automatic=make_mcmc_grouping(edges,lines);grouping=&automatic; }
    // Replay complete boundary increments, including post-absorption bond shifts.
    // Environment construction and local update plans consume this ordered graph;
    // never reconstruct it using absorption alone.
    target_edges = edges;
    std::size_t min_line_idx = 0;
    std::size_t max_line_idx = (lines.size() > 0) ? (lines.size()-1) : 0;
    std::size_t target_line_idx = static_cast<std::size_t>(target_line);
    if( target_line_idx > 0 ) {
      for(std::size_t k=min_line_idx+1; k < target_line_idx; k++) {
	bmps_increment_line_for_edges(lines[k-1],lines[k],target_edges,&(*grouping)[k]);
      }
    }
    for(std::size_t k=max_line_idx; k > target_line_idx+1; k--) {
      bmps_increment_line_for_edges(lines[k],lines[k-1],target_edges,&(*grouping)[k-1]);
    }
  }

  /**
   * @brief Construct the initial lower and upper boundary MPS caches.
   * @tparam TenT TCAPI tensor type.
   * @tparam IntT Signed integer type for site labels, states and indices.
   * @param[in,out] ctx TCAPI tensor backend context.
   * @param[in] edges Ordered undirected graph edges; virtual axes follow incident-edge order.
   * @param[in] V Attached local tensors in site_idx order.
   * @param[in] site_idx Site labels in the same local order as the tensor/state arrays.
   * @param[in] site_to_mpi_rank Global site-to-owner-rank map; original ownership is unchanged.
   * @param[in] lines Ordered nonempty lines, from lower to upper boundary.
   * @param[out] T_upper Boundary tensors indexed [line][local tensor]; upper region absorbed onto that line.
   * @param[out] site_idx_upper Local site labels for each T_upper[line].
   * @param[in,out] E_upper BP messages indexed [line][local message], initialized where BP is run.
   * @param[in,out] edge_idx_upper Local BP edge indices per line, relative to that boundary induced graph.
   * @param[out] T_lower Boundary tensors indexed [line][local tensor]; lower region absorbed onto that line.
   * @param[out] site_idx_lower Local site labels for each T_lower[line].
   * @param[in,out] E_lower BP messages indexed [line][local message], initialized where BP is run.
   * @param[in,out] edge_idx_lower Local BP edge indices per line, relative to that boundary induced graph.
   * @param[in] max_bp_steps Maximum BP iterations per boundary truncation.
   * @param[in] bp_tolerance BP residual stopping tolerance.
   * @param[in] max_bond_dim Maximum retained boundary bond dimension.
   * @param[in] sv_min Singular-value threshold passed to boundary truncation.
   * @param[in] tg_err Target truncation error passed to boundary truncation.
   * @param[in] comm MPI communicator; all ranks must participate in the same call order.
   * @param[in] grouping Optional fixed three-line grouping. Null selects automatic grouping from the original graph.
   *
   * T_lower[k] represents lines 0 through k; T_upper[k] represents lines k through
   * the last line. A target k consumes lower[k-1] and upper[k+1] when present.
   * Only needed boundaries are built: lower at the last line and upper at line zero
   * are unused. BP messages are meaningful only where absorption/truncation runs.
   * @pre Lines and grouping form a valid sampling plan; V and owner metadata agree.
   * Collective tensor operations use the saved groups in both directions.
   * @see bmps_increment_line
   * @see bmps_bp_truncation
   */
  template <typename TenT, typename IntT>
  void mcmc_init_bmps(context_handle_t<TenT> & ctx,
		      const std::vector<std::pair<IntT,IntT>> & edges,
		      const std::vector<TenT> & V,
		      const std::vector<IntT> & site_idx,
		      const std::map<IntT,IntT> & site_to_mpi_rank,
		      const std::vector<std::vector<IntT>> & lines,
		      std::vector<std::vector<TenT>> & T_upper,
		      std::vector<std::vector<IntT>> & site_idx_upper,
		      std::vector<std::vector<TenT>> & E_upper,
		      std::vector<std::vector<IntT>> & edge_idx_upper,
		      std::vector<std::vector<TenT>> & T_lower,
		      std::vector<std::vector<IntT>> & site_idx_lower,
		      std::vector<std::vector<TenT>> & E_lower,
		      std::vector<std::vector<IntT>> & edge_idx_lower,
		      IntT max_bp_steps,
		      real_t<TenT> bp_tolerance,
		      bond_dim_t<TenT> max_bond_dim,
		      real_t<TenT> sv_min,
		      real_t<TenT> tg_err,
		      MPI_Comm comm,
       const MCMCGrouping<IntT>* grouping=nullptr) {
    MCMCGrouping<IntT> automatic;
    if(!grouping) { automatic=make_mcmc_grouping(edges,lines);grouping=&automatic; }


    using RealT = typename tcapi::tensor_traits<TenT>::real_t;
    using BondDimT = typename tcapi::tensor_traits<TenT>::bond_dim_t;

    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(comm,&mpi_size);

    std::size_t min_line_idx = 0;
    std::size_t max_line_idx = (lines.size() > 1) ? (lines.size()-1) : 0;

    // bmps from lower edges
    T_lower.resize(lines.size());
    site_idx_lower.resize(lines.size());
    E_lower.resize(lines.size());
    edge_idx_lower.resize(lines.size());
    // No target line consumes a lower boundary at the uppermost line.
    // Clear the unused slot as well when output containers are reused.
    if(!lines.empty()) {
      T_lower.back().clear();
      site_idx_lower.back().clear();
      E_lower.back().clear();
      edge_idx_lower.back().clear();
    }
    std::vector<std::pair<IntT,IntT>> edges_lower(edges);
    for(std::size_t target_line_idx=min_line_idx;
	target_line_idx < max_line_idx; target_line_idx++) {
      T_lower[target_line_idx].resize(0);
      site_idx_lower[target_line_idx].resize(0);
      if( target_line_idx > min_line_idx ) {
	for(std::size_t adrs=0; adrs < lines[target_line_idx-1].size(); adrs++) {
	  auto site_this = lines[target_line_idx-1][adrs];
	  auto mpi_rank_this = site_to_mpi_rank.at(site_this);
	  if( mpi_rank_this == mpi_rank ) {
	    auto it_site_adrs_this = std::find(site_idx_lower[target_line_idx-1].begin(),
					       site_idx_lower[target_line_idx-1].end(),
					       site_this);
	    auto site_adrs_this = std::distance(site_idx_lower[target_line_idx-1].begin(),
						it_site_adrs_this);
	    T_lower[target_line_idx].push_back(
	      tcapi::copy(ctx,T_lower[target_line_idx-1][site_adrs_this]));
	    site_idx_lower[target_line_idx].push_back(site_this);
	  }
	}
      }
      for(std::size_t adrs=0; adrs < lines[target_line_idx].size(); adrs++) {
	auto site_this = lines[target_line_idx][adrs];
	auto mpi_rank_this = site_to_mpi_rank.at(site_this);
	if( mpi_rank_this == mpi_rank ) {
	  auto it_site_adrs_this = std::find(site_idx.begin(),site_idx.end(),site_this);
	  auto site_adrs_this = std::distance(site_idx.begin(),it_site_adrs_this);
	  T_lower[target_line_idx].push_back(tcapi::copy(ctx,V[site_adrs_this]));
	  site_idx_lower[target_line_idx].push_back(site_this);
	}
      }
      if( target_line_idx > min_line_idx ) {
	std::vector<BondDimT> res_bond_dim;
	std::vector<RealT> res_trunc_err;
	auto site_to_mpi_rank_temporary = site_to_mpi_rank;
	bmps_increment_line(ctx,
			    lines[target_line_idx-1],
			    lines[target_line_idx],
			    edges_lower,
			    T_lower[target_line_idx],
			    site_idx_lower[target_line_idx],
			    site_to_mpi_rank_temporary,
			    comm, &(*grouping)[target_line_idx]);
	// perform bp and truncation
	bmps_bp_truncation(ctx,
			   edges_lower,lines[target_line_idx],
			   T_lower[target_line_idx],site_idx_lower[target_line_idx],
			   site_to_mpi_rank_temporary,
			   E_lower[target_line_idx],edge_idx_lower[target_line_idx],
			   max_bp_steps,bp_tolerance,
			   max_bond_dim,sv_min,tg_err,
			   res_bond_dim,res_trunc_err,
			   comm,true);
      }
    }

#ifdef MCMC_DEBUG
    if( mpi_rank == 0 ) {
      std::cout << " End construction of lower bmps in mcmc_init_bmps" << std::endl;
    }
#endif

    // bmps from upper edges
    T_upper.resize(lines.size());
    site_idx_upper.resize(lines.size());
    E_upper.resize(lines.size());
    edge_idx_upper.resize(lines.size());
    std::vector<std::pair<IntT,IntT>> edges_upper(edges);
    for(std::size_t target_line_idx=max_line_idx;
	target_line_idx > min_line_idx; target_line_idx--) {
      T_upper[target_line_idx-min_line_idx].resize(0);
      site_idx_upper[target_line_idx-min_line_idx].resize(0);
      if( target_line_idx < max_line_idx ) {
	for(std::size_t adrs=0; adrs < lines[target_line_idx+1].size(); adrs++) {
	  auto site_this = lines[target_line_idx+1][adrs];
	  auto mpi_rank_this = site_to_mpi_rank.at(site_this);
	  if( mpi_rank_this == mpi_rank ) {
	    auto it_site_adrs_this = std::find(site_idx_upper[target_line_idx+1].begin(),
					       site_idx_upper[target_line_idx+1].end(),
					       site_this);
	    auto site_adrs_this = std::distance(site_idx_upper[target_line_idx+1].begin(),
						it_site_adrs_this);
	    T_upper[target_line_idx].push_back(
	      tcapi::copy(ctx,T_upper[target_line_idx+1][site_adrs_this]));
	    site_idx_upper[target_line_idx].push_back(site_this);
	  }
	}
      }
      for(std::size_t adrs=0; adrs < lines[target_line_idx].size(); adrs++) {
	auto site_this = lines[target_line_idx][adrs];
	auto mpi_rank_this = site_to_mpi_rank.at(site_this);
	if( mpi_rank_this == mpi_rank ) {
	  auto it_site_adrs_this = std::find(site_idx.begin(),site_idx.end(),site_this);
	  auto site_adrs_this = std::distance(site_idx.begin(),it_site_adrs_this);
	  T_upper[target_line_idx].push_back(tcapi::copy(ctx,V[site_adrs_this]));
	  site_idx_upper[target_line_idx].push_back(site_this);
	}
      }
      if( target_line_idx < max_line_idx ) {
	std::vector<BondDimT> res_bond_dim;
	std::vector<RealT> res_trunc_err;
	auto site_to_mpi_rank_temporary = site_to_mpi_rank;
	bmps_increment_line(ctx,
			    lines[target_line_idx+1],
			    lines[target_line_idx],
			    edges_upper,
			    T_upper[target_line_idx],
			    site_idx_upper[target_line_idx],
			    site_to_mpi_rank_temporary,
			    comm, &(*grouping)[target_line_idx]);
	// perform bp and truncation
	bmps_bp_truncation(ctx,edges_upper,lines[target_line_idx],
			   T_upper[target_line_idx],site_idx_upper[target_line_idx],
			   site_to_mpi_rank_temporary,
			   E_upper[target_line_idx],edge_idx_upper[target_line_idx],
			   max_bp_steps,bp_tolerance,
			   max_bond_dim,sv_min,tg_err,
			   res_bond_dim,res_trunc_err,
			   comm,true);
      }
    }
#ifdef MCMC_DEBUG
    if( mpi_rank == 0 ) {
      std::cout << " End construction of upper bmps in mcmc_init_bmps" << std::endl;
    }
#endif
  }

  /**
   * @brief Build proper prefix or suffix environments from saved groups.
   * @tparam TenT TCAPI tensor type.
   * @tparam IntT Signed integer type for site labels, states and indices.
   * @param[in,out] ctx TCAPI tensor backend context.
   * @param[in] site_to_mpi_rank Global site-to-owner-rank map; original ownership is unchanged.
   * @param[in] edges Ordered undirected graph edges; virtual axes follow incident-edge order.
   * @param[in] line_upper Upper neighboring boundary sites; empty at the upper edge.
   * @param[in] T_upper Rank-local upper boundary tensors.
   * @param[in] site_idx_upper Labels corresponding to local T_upper entries.
   * @param[in] line_target Target-line sites in left-to-right order.
   * @param[in] T_target Rank-local attached tensors of the target line.
   * @param[in] site_idx_target Labels corresponding to local T_target entries.
   * @param[in] line_lower Lower neighboring boundary sites; empty at the lower edge.
   * @param[in] T_lower Rank-local lower boundary tensors.
   * @param[in] site_idx_lower Labels corresponding to local T_lower entries.
   * @param[out] Env Local cumulative environment tensors, cleared before construction.
   * @param[out] site_idx_env Target-site labels identifying the stored environments, in local construction order.
   * @param[in] direction 0 builds from left to right; nonzero builds from right to left.
   * @param[in] site_line_sets Fixed groups, one per target site, including adjacent boundary members.
   * @param[in] comm MPI communicator; all ranks must participate in the same call order.
   *
   * For N target sites, build N-1 environments globally; N <= 1 returns empty
   * outputs. The full-line scalar contraction is not stored. Each cumulative
   * environment belongs to the owner of its last included target site and is keyed
   * by that target label. Groups are followed as supplied, never recomputed here.
   * @pre Tensors, edges and groups describe the same reduced three-line network.
   * @see mcmc_update_line
   */
  template <typename TenT, typename IntT>
  void mcmc_init_environment(context_handle_t<TenT> & ctx,
			     const std::map<IntT,IntT> & site_to_mpi_rank,
			     const std::vector<std::pair<IntT,IntT>> & edges,
			     const std::vector<IntT> & line_upper,
			     const std::vector<TenT> & T_upper,
			     const std::vector<IntT> & site_idx_upper,
			     const std::vector<IntT> & line_target,
			     const std::vector<TenT> & T_target,
			     const std::vector<IntT> & site_idx_target,
			     const std::vector<IntT> & line_lower,
			     const std::vector<TenT> & T_lower,
			     const std::vector<IntT> & site_idx_lower,
			     std::vector<TenT> & Env,
			     std::vector<IntT> & site_idx_env,
			     IntT direction,
			     const std::vector<std::vector<IntT>> & site_line_sets,
			     MPI_Comm comm) {

    // Only proper prefixes/suffixes are needed by the next site's update.
    // The full-line contraction is never used as an opposite environment.
    Env.clear();
    site_idx_env.clear();
    if(line_target.size() <= 1) return;

    // Follow the saved groups; do not select neighbors again in this direction.
    // Store each cumulative environment on the owner of its last target site.

    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(comm,&mpi_size);

    std::size_t min_adrs_site_line = 0;
    std::size_t max_adrs_site_line = ( line_target.size() > 1 ) ? (line_target.size()-1) : 0;

    auto all_site_label = GetSiteIndexFromBond(edges);
    IntT nenv_site_label = (all_site_label.size() > 0) ?
      (all_site_label[all_site_label.size()-1]+1) : -2;
    IntT tenv_site_label = nenv_site_label+1;


#ifdef MCMC_DEBUG
    if( mpi_rank == 0 ) {
      std::cout << " mcmc_init_environment: input edges =:";
      for(const auto & [a,b] : edges) {
	std::cout << " (" << a << "," << b << ")";
      }
      std::cout << std::endl;
    }
#endif

    // we first define the sets of sites consisting of the environments
#ifdef MCMC_DEBUG
    if( mpi_rank == 0 ) {
      std::cout << " mcmc_init_environment: site blocks =";
      for(const auto & sites : site_line_sets) {
	std::cout << " (";
	for(const auto & site : sites) {
	  std::cout << " " << site;
	}
	std::cout << ")";
      }
      std::cout << std::endl;
    }
#endif

    // construct the left environment tensors using environment tensors
    auto edges_temporary = edges;
    std::vector<TenT> T;
    std::vector<IntT> site_idx;
    T.resize(0);
    site_idx.resize(0);
    for(std::size_t adrs_site_line = 0;
	adrs_site_line < site_idx_lower.size();
	adrs_site_line++) {
      auto site = site_idx_lower[adrs_site_line];
      auto rank = site_to_mpi_rank.at(site);
      if( rank == mpi_rank ) {
	T.push_back(tcapi::copy(ctx,T_lower[adrs_site_line]));
	site_idx.push_back(site);
      }
    }
    for(std::size_t adrs_site_line = 0;
	adrs_site_line < site_idx_target.size();
	adrs_site_line++) {
      auto site = site_idx_target[adrs_site_line];
      auto rank = site_to_mpi_rank.at(site);
      if( rank == mpi_rank ) {
	T.push_back(tcapi::copy(ctx,T_target[adrs_site_line]));
	site_idx.push_back(site);
      }
    }
    for(std::size_t adrs_site_line = 0;
	adrs_site_line < site_idx_upper.size();
	adrs_site_line++) {
      auto site = site_idx_upper[adrs_site_line];
      auto rank = site_to_mpi_rank.at(site);
      if( rank == mpi_rank ) {
	T.push_back(tcapi::copy(ctx,T_upper[adrs_site_line]));
	site_idx.push_back(site);
      }
    }
#ifdef MCMC_DEBUG
    if( mpi_rank == 0 ) {
      std::cout << " mcmc_init_environment: start the construction for the environment tensors, target =";
      for(const auto & site : line_target) {
	std::cout << " " << site << "(rank " << site_to_mpi_rank.at(site) << ")";
      }
      std::cout << std::endl;
      std::cout << " mcmc_init_environment: upper lines =";
      for(const auto & site : line_upper) {
	std::cout << " " << site << "(rank " << site_to_mpi_rank.at(site) << ")";
      }
      std::cout << std::endl;
      std::cout << " mcmc_init_environment: lower lines =";
      for(const auto & site : line_lower) {
	std::cout << " " << site << "(rank " << site_to_mpi_rank.at(site) << ")";
      }
      std::cout << std::endl;
    }
    /*
    for(int rank=0; rank < mpi_size; rank++) {
      if( rank == mpi_rank ) {
	std::cout << " mcmc_init_environment: site_idx at rank "
		  << mpi_rank << " =";
	for(std::size_t adrs=0; adrs < site_idx.size(); adrs++) {
	  std::cout << " " << site_idx[adrs];
	}
      }
      MPI_Barrier(comm);
    }
    std::cout << std::endl;
    */
#endif
    std::map<IntT,IntT> site_to_mpi_rank_temp;
    for(const auto & site : line_lower) {
      site_to_mpi_rank_temp[site] = site_to_mpi_rank.at(site);
    }
    for(const auto & site : line_target) {
      site_to_mpi_rank_temp[site] = site_to_mpi_rank.at(site);
    }
    for(const auto & site : line_upper) {
      site_to_mpi_rank_temp[site] = site_to_mpi_rank.at(site);
    }
    int pos_begin = 0;
    int pos_end   = max_adrs_site_line - min_adrs_site_line;
    int adrs_site_line_begin = (direction==0) ?
      min_adrs_site_line : max_adrs_site_line;
    int diff_pos = (direction==0) ? 1 : -1;
    // Exclude the final site in either direction: build exactly N-1 environments.
    for(int pos = pos_begin; pos < pos_end; pos++) {
      auto adrs_site_line = adrs_site_line_begin + diff_pos * pos;
      auto sites_contract = site_line_sets[adrs_site_line];
      if( adrs_site_line != adrs_site_line_begin ) {
	sites_contract.push_back(nenv_site_label);
      }
      auto edges_contract = extract_induced_edges(edges_temporary,
						  sites_contract);
      std::map<IntT,IntT> site_priority;
      std::map<IntT,IntT> label_map;
      site_priority[line_target[adrs_site_line]] = 1;
      label_map[line_target[adrs_site_line]] = tenv_site_label;
#ifdef MCMC_DEBUG
      int mpi_rank_target = site_to_mpi_rank.at(line_target[adrs_site_line]);
      if( mpi_rank == mpi_rank_target ) {
	std::cout << " mcmc_init_environment: start the contraction for sites";
	for(auto const & x : sites_contract) {
	  std::cout << " " << x << " (rank " << site_to_mpi_rank_temp.at(x) << ")";
	}
	std::cout << ", edges =";
	for(const auto & [u,v] : edges_temporary) {
	  std::cout << " (" << u << "," << v << ")";
	}
	std::cout << std::endl;
	std::cout << " mcmc_init_environment: edges will be contracted =";
	for(const auto & [u,v] : edges_contract) {
	  std::cout << " (" << u << "," << v << ")";
	}
	std::cout << std::endl;
      }
#endif
      if(edges_contract.empty()) {
        if(sites_contract.size()!=1)
          throw std::invalid_argument("environment group has no contraction path");
        graph_rename_site(edges_temporary,site_idx,site_to_mpi_rank_temp,
                          sites_contract.front(),tenv_site_label);
      }
      graph_tensor_network_contraction(ctx,edges_temporary,
				       T,site_idx,
				       site_to_mpi_rank_temp,
				       edges_contract,
				       site_priority,
				       label_map,comm);
#ifdef MCMC_DEBUG
      if( mpi_rank == mpi_rank_target ) {
	std::cout << " mcmc_init_environment: end the contractions for local sites";
	for(auto const & x : site_line_sets[adrs_site_line]) {
	  if( x == line_target[adrs_site_line] ) {
	    std::cout << " " << x
		      << "->" << tenv_site_label
		      << " (rank " << site_to_mpi_rank_temp.at(tenv_site_label)
		      << ")";

	  } else {
	    std::cout << " " << x
		      << " (rank none)";
	  }
	}
	std::cout << ", edges = ";
	for(const auto & [u,v] : edges_temporary) {
	  std::cout << " (" << u << "," << v << ")";
	}
	std::cout << std::endl;
      }
#endif
      for(std::size_t adrs_edge=0; adrs_edge < edges_temporary.size(); adrs_edge++) {
	if( edges_temporary[adrs_edge].first == tenv_site_label ) {
	  auto site_f = nenv_site_label;
	  auto site_s = edges_temporary[adrs_edge].second;
	  edges_temporary[adrs_edge] = make_edge(site_f,site_s);
	} else if ( edges_temporary[adrs_edge].second == tenv_site_label ) {
	  auto site_f = edges_temporary[adrs_edge].first;
	  auto site_s = nenv_site_label;
	  edges_temporary[adrs_edge] = make_edge(site_f,site_s);
	}
      }
      if( site_to_mpi_rank_temp.at(tenv_site_label) == mpi_rank ) {
	auto it_adrs_env = std::find(site_idx.begin(),
				     site_idx.end(),
				     tenv_site_label);
	auto adrs_env = std::distance(site_idx.begin(),
				      it_adrs_env);
	site_idx[adrs_env] = nenv_site_label;
	Env.push_back(tcapi::copy(ctx,T[adrs_env]));
	site_idx_env.push_back(line_target[adrs_site_line]);
      }
      auto it_prev_nenv_rank = site_to_mpi_rank_temp.find(nenv_site_label);
      if( it_prev_nenv_rank != site_to_mpi_rank_temp.end() ) {
	site_to_mpi_rank_temp.erase(it_prev_nenv_rank);
      }
      auto node = site_to_mpi_rank_temp.extract(tenv_site_label);
      if (!node.empty()) {
	node.key() = nenv_site_label;
	site_to_mpi_rank_temp.insert(std::move(node));
      }
    }
  }

  /**
   * @brief Adjust local tensor axes to a graph with grouped left/right environments.
   * @tparam TenT TCAPI tensor type.
   * @tparam IntT Signed integer type for site labels, states and indices.
   * @param[in,out] ctx TCAPI tensor backend context.
   * @param[in] line Ordered nonempty target line.
   * @param[in] edges Ordered undirected graph edges; virtual axes follow incident-edge order.
   * @param[in] site_env Saved groups indexed by target-line position.
   * @param[in,out] T Supplied local tensors whose axes must follow environment edge fusion.
   * @param[in] site_idx Site labels in the same local order as the tensor/state arrays.
   * @param[in] site_to_mpi_rank Owners of supplied tensor sites; absent graph sites are assigned internal owner -1.
   * @param[in] left_env_adrs Last group included in the left environment; -1 skips the left side.
   * @param[in] right_env_adrs First group included in the right environment; line.size() skips the right side.
   * @param[out] edges_result Ordered graph after the specified environment contractions.
   * @param[in] comm MPI communicator; all ranks must participate in the same call order.
   *
   * Uses the internal partial-contraction path for absent tensors and axis fusion
   * on supplied affected tensors. Local label and owner metadata are copied internally.
   * This is preparation for a local update, not construction of all environment values.
   * @pre IntT is signed, group bounds are valid, and supplied tensors/owners meet
   * the partial-contraction contract.
   * @see mcmc_update_line
   */
  template <typename TenT, typename IntT>
  void mcmc_prep_local_target_tensors(context_handle_t<TenT> & ctx,
				      const std::vector<IntT> & line,
				      const std::vector<std::pair<IntT,IntT>> & edges,
				      const std::vector<std::vector<IntT>> & site_env,
				      std::vector<TenT> & T,
				      const std::vector<IntT> & site_idx,
				      const std::map<IntT,IntT> site_to_mpi_rank,
				      IntT left_env_adrs,
				      IntT right_env_adrs,
				      std::vector<std::pair<IntT,IntT>> & edges_result,
				      MPI_Comm comm) {
    // assume that T is only defined at one of site env position
    // left_env_adrs < 0 indicates no Lenv position
    // right_env_adrs > site_env.size()-1 indicates no Renv position

    edges_result = edges;
    auto site_idx_temp = site_idx;
    auto all_sites = GetSiteIndexFromBond(edges_result);
    std::map<IntT,IntT> site_to_mpi_rank_temp;
    for(const auto & site : all_sites) {
      auto it = site_to_mpi_rank.find(site);
      if( it == site_to_mpi_rank.end() ) {
	site_to_mpi_rank_temp[site] = -1;
      } else {
	site_to_mpi_rank_temp[site] = site_to_mpi_rank.at(site);
      }
    }
    for(IntT pos=0; pos <= left_env_adrs; pos++) {
      std::vector<IntT> sites = site_env[pos];
      if( pos > 0 ) {
	sites.push_back(line[pos-1]);
      }
      auto edges_contract = extract_induced_edges(edges_result,sites);
      std::map<IntT,IntT> site_priority;
      std::map<IntT,IntT> label_map;
      site_priority[line[pos]] = 1;
      detail::graph_tensor_network_partial_contraction(ctx,edges_result,
				       T,site_idx_temp,
				       site_to_mpi_rank_temp,
				       edges_contract,
				       site_priority,label_map,comm);
    }

    IntT num_line_end = static_cast<IntT>(line.size()-1);

    for(IntT pos=num_line_end; pos >= right_env_adrs; pos--) {
      std::vector<IntT> sites = site_env[pos];
      if( pos < num_line_end ) {
	sites.push_back(line[pos+1]);
      }
      auto edges_contract = extract_induced_edges(edges_result,sites);
      std::map<IntT,IntT> site_priority;
      std::map<IntT,IntT> label_map;
      site_priority[line[pos]] = 1;
      detail::graph_tensor_network_partial_contraction(ctx,edges_result,
				       T,site_idx_temp,
				       site_to_mpi_rank_temp,
				       edges_contract,
				       site_priority,
				       label_map,comm);
    }
  }


  /**
   * @brief Update one-/two-site windows on a line using boundary and side environments.
   * @tparam TenT TCAPI tensor type.
   * @tparam IntT Signed integer type for site labels, states and indices.
   * @tparam RNG Uniform random bit generator type.
   * @param[in,out] ctx TCAPI tensor backend context.
   * @param[in] edges Ordered undirected graph edges; virtual axes follow incident-edge order.
   * @param[in] W Original rank-local tensors, with the physical axis last; never modified.
   * @param[in] site_idx Site labels in the same local order as the tensor/state arrays.
   * @param[in] site_to_mpi_rank Global site-to-owner-rank map; original ownership is unchanged.
   * @param[in] line_upper Upper neighboring boundary sites; empty at the upper edge.
   * @param[in] T_upper Rank-local upper boundary tensors.
   * @param[in] site_idx_upper Labels corresponding to local T_upper entries.
   * @param[in] line_target Target-line sites in left-to-right order.
   * @param[in,out] T_target Attached target-line tensors, replaced on accepted proposals.
   * @param[in] site_idx_target Labels corresponding to local T_target entries.
   * @param[in] line_lower Lower neighboring boundary sites; empty at the lower edge.
   * @param[in] T_lower Rank-local lower boundary tensors.
   * @param[in] site_idx_lower Labels corresponding to local T_lower entries.
   * @param[in] site_line_sets Saved groups corresponding to line_target; never regrouped here.
   * @param[in,out] site_states Current physical states in rank-local site_idx order.
   * @param[in] comm MPI communicator; all ranks must participate in the same call order.
   * @param[in,out] rng Rank-local generator; method 2 draws on the leading target owner; methods 0/1 on site owners.
   * @param[in] sampling_method 0 for Metropolis; 1 for Barker; 2 for full-state heat bath.
   * @param[in] update_sites Window size 1 or 2; size 2 requires method 2.
   * @param[in] sampling_direction 0 visits left to right; nonzero visits right to left.
   *
   * Method 2 delegates full-state windows to mcmc_update_line_heat_bath.
   * Methods 0/1 first build opposite-side environments, then advance the current-side environment
   * using each accepted or retained state. Propose a different state symmetrically
   * and evaluate current/proposed amplitudes through the same contraction path.
   * Physical dimension one is fixed and consumes no proposal or acceptance random
   * numbers. Original W and upper/lower boundary inputs are not modified.
   * @pre The target line is nonempty; tensors, states and saved groups are consistent.
   * All ranks enter collectively. Boundary truncation can make probabilities approximate.
   * @throws std::invalid_argument If the sampling method or a delegated operation rejects input.
   * @see mcmc_init_environment
   * @see mcmc_pair_metropolis_probability
   * @see mcmc_pair_heat_bath_probability
   */
  template <typename TenT, typename IntT, typename RNG>
  void mcmc_update_line(context_handle_t<TenT> & ctx,
			const std::vector<std::pair<IntT,IntT>> & edges,
			const std::vector<TenT> & W,
			const std::vector<IntT> & site_idx,
			const std::map<IntT,IntT> & site_to_mpi_rank,
			const std::vector<IntT> & line_upper,
			const std::vector<TenT> & T_upper,
			const std::vector<IntT> & site_idx_upper,
			const std::vector<IntT> & line_target,
			std::vector<TenT> & T_target,
			const std::vector<IntT> & site_idx_target,
			const std::vector<IntT> & line_lower,
			const std::vector<TenT> & T_lower,
			const std::vector<IntT> & site_idx_lower,
			const std::vector<std::vector<IntT>> & site_line_sets,
			std::vector<IntT> & site_states,
			MPI_Comm comm,
			RNG & rng,
			int sampling_method = 2,
			int sampling_direction = 0, int update_sites = 1,
                        const std::map<std::vector<int>,MPI_Comm>* region_communicators=nullptr,
                        const mcmc_edge_line_plan<IntT>* edge_plan=nullptr) {
    if(update_sites!=1 && update_sites!=2) throw std::invalid_argument("update_sites must be 1 or 2");
    if(sampling_method==2) {
      mcmc_update_line_heat_bath(ctx,edges,W,site_idx,site_to_mpi_rank,
        T_upper,site_idx_upper,T_target,site_idx_target,T_lower,site_idx_lower,
        line_target,site_line_sets,site_states,comm,rng,sampling_direction,update_sites,nullptr,region_communicators,edge_plan);
      return;
    }
    if(update_sites!=1) throw std::invalid_argument("two-site updates require full-state heat bath (method 2)");

    using ElemT = typename tcapi::tensor_traits<TenT>::elem_t;
    using RealT = typename tcapi::tensor_traits<TenT>::real_t;
    using BondLabelT = typename tcapi::tensor_traits<TenT>::bond_label_t;
    using ShapeT = typename tcapi::tensor_traits<TenT>::shape_t;
    using CoorsT = typename tcapi::tensor_traits<TenT>::elem_coors_t;

    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    std::uniform_real_distribution<RealT> dist_real(0.0,1.0);
#ifdef MCMC_DEBUG
    if( mpi_rank == 0 ) {
      std::cout << " mcmc_update_line: start ================================================ "
		<< std::endl;
    }
#endif
    validate_mcmc_sampling_method(sampling_method);
    // Both methods share the symmetric single-site proposal and contractions.
    {

      TenT Nenv;
      IntT site_label_nenv;
      int mpi_rank_nenv;
      std::vector<TenT> T_tenv;
      std::vector<IntT> site_idx_tenv;
      IntT tenv_direction = (sampling_direction == 0) ? 1 : 0;
      mcmc_init_environment(ctx,site_to_mpi_rank,edges,
                            line_upper,T_upper,site_idx_upper,
                            line_target,T_target,site_idx_target,
                            line_lower,T_lower,site_idx_lower,
                            T_tenv,site_idx_tenv,tenv_direction,
                            site_line_sets,comm);
#ifdef MCMC_DEBUG
      std::vector<int> send_sync_message(mpi_size,0);
      std::vector<int> recv_sync_message(mpi_size,0);
      send_sync_message[mpi_rank] = 1;
      MPI_Allreduce(send_sync_message.data(),
		    recv_sync_message.data(),
		    mpi_size,MPI_INT,MPI_SUM,comm);
      if ( mpi_rank == 0 ) {
	std::cout << " mcmc_update_line: ---------------------------- end mcmc_init_environment "
		  << std::endl;
      }
      MPI_Barrier(comm);
#endif

      // we first construct the set of edges for the local contractions
      auto all_site_label = GetSiteIndexFromBond(edges);
      IntT lenv_site_label = (all_site_label.size() > 0) ?
	(all_site_label[all_site_label.size()-1]+1) : -2;
      IntT renv_site_label = lenv_site_label+1;

      // Start the contraction from the left/right to the right/left
      std::size_t min_local_adrs = 0;
      std::size_t max_local_adrs = (site_line_sets.size() > 0) ?
	(site_line_sets.size()-1) : 0;
      int local_adrs_begin;
      int local_adrs_end;
      int local_adrs_diff;
      if( sampling_direction == 0 ) {
	local_adrs_begin = static_cast<int>(min_local_adrs);
	local_adrs_end   = static_cast<int>(max_local_adrs);
	local_adrs_diff  = 1;
      } else {
	local_adrs_begin = static_cast<int>(max_local_adrs);
	local_adrs_end   = static_cast<int>(min_local_adrs);
	local_adrs_diff  = -1;
      }
      for(int local_point = 0;
	  local_point < static_cast<int>(site_line_sets.size());
	  local_point++) {
	int local_adrs = local_adrs_begin + local_adrs_diff * local_point;
#ifdef MCMC_DEBUG
	if( mpi_rank == 0 ) {
	  std::cout << " mcmc_update_line: ---------- start local update at site "
		    << line_target[local_adrs]
		    << " ---------- " << std::endl;
	}
#endif
	auto local_sites = site_line_sets[local_adrs];
	std::vector<TenT> T_local;
	std::vector<IntT> site_idx_local;
	std::map<IntT,IntT> site_to_mpi_rank_local;
	for(std::size_t site_adrs=0;
	    site_adrs < local_sites.size();
	    site_adrs++) {
	  auto site_label_this = local_sites[site_adrs];
	  auto mpi_rank_this = site_to_mpi_rank.at(site_label_this);
	  site_to_mpi_rank_local[site_label_this] = mpi_rank_this;
	  if( mpi_rank == mpi_rank_this ) {
	    auto it_adrs_upper = std::find(site_idx_upper.begin(),site_idx_upper.end(),
					   site_label_this);
	    auto it_adrs_target = std::find(site_idx_target.begin(),site_idx_target.end(),
					    site_label_this);
	    auto it_adrs_lower = std::find(site_idx_lower.begin(),site_idx_lower.end(),
					   site_label_this);
	    if( it_adrs_lower != site_idx_lower.end() ) {
	      auto adrs_lower = std::distance(site_idx_lower.begin(),
					      it_adrs_lower);
	      T_local.push_back(tcapi::copy(ctx,T_lower[adrs_lower]));
	      site_idx_local.push_back(site_label_this);
	    } else if ( it_adrs_upper != site_idx_upper.end() ) {
	      auto adrs_upper = std::distance(site_idx_upper.begin(),
					      it_adrs_upper);
		T_local.push_back(tcapi::copy(ctx,T_upper[adrs_upper]));
		site_idx_local.push_back(site_label_this);
	      } else if ( it_adrs_target != site_idx_target.end() ) {
	      auto adrs_target = std::distance(site_idx_target.begin(),
					       it_adrs_target);
	      T_local.push_back(tcapi::copy(ctx,T_target[adrs_target]));
	      site_idx_local.push_back(site_label_this);
	    }
	  }
	} // end for(std::size_t site_adrs)
#ifdef MCMC_DEBUG
	for(int rank=0; rank < mpi_size; rank++) {
	  if( mpi_rank == rank ) {
	    std::cout << " mcmc_update_line: end preparation of T_local (size="
		      << T_local.size() << ") at rank " << mpi_rank;
	    for(const auto & site : site_idx_local) {
	      std::cout << ", site " << site
			<< " (mpi " << site_to_mpi_rank_local.at(site) << ")";
	    }
	    std::cout << std::endl;
	  }
	  MPI_Barrier(comm);
	}
	MPI_Barrier(comm);
#endif
	std::vector<TenT> T_update;
        T_update.reserve(T_local.size());
        for(const auto& t : T_local) T_update.push_back(tcapi::copy(ctx,t));
	std::vector<IntT> site_idx_update(site_idx_local);
	std::map<IntT,IntT> site_to_mpi_rank_update(site_to_mpi_rank_local);
	auto mpi_rank_target = site_to_mpi_rank.at(line_target[local_adrs]);
	int state_update;
        bool fixed_site = false; // Used only by the target site's owning rank.
	TenT V_update;
	if( mpi_rank == mpi_rank_target ) {
	  auto it_site_local_adrs = std::find(site_idx_local.begin(),
					      site_idx_local.end(),
					      line_target[local_adrs]);
	  auto site_local_adrs = std::distance(site_idx_local.begin(),
					       it_site_local_adrs);
	  auto it_site_adrs = std::find(site_idx.begin(),
					site_idx.end(),
					line_target[local_adrs]);
	  auto site_adrs = std::distance(site_idx.begin(),
					 it_site_adrs);
	  auto shape_w = tcapi::shape(ctx,W[site_adrs]);
	  auto order_w = tcapi::order(ctx,W[site_adrs]);
          if(order_w == 0 || shape_w[order_w-1] < 1)
            throw std::invalid_argument("sampling requires a nonempty physical index");
          int d = static_cast<int>(shape_w[order_w-1]);
          fixed_site = (d == 1);
          state_update = site_states[site_adrs];
          if(fixed_site) {
            if(state_update != 0)
              throw std::invalid_argument("a dimension-one site's state must be zero");
          } else {
            std::uniform_int_distribution<int> dist(0,d-2);
            state_update = dist(rng);
            if(state_update >= site_states[site_adrs]) ++state_update;
          }
          // Keep the existing attached tensor for a fixed site. The common
          // contraction path below must still advance the sweep environment.
          if(!fixed_site) {
	  TenT S;
	  state_tensor(ctx,state_update,d,S);
	  auto order_s = tcapi::order(ctx,S);
	  auto order_t = tcapi::order(ctx,T_local[site_local_adrs]);
	  tcapi::List<BondLabelT> label_w(order_w);
	  tcapi::List<BondLabelT> label_s(order_s);
	  tcapi::List<BondLabelT> label_t(order_t);
	  std::iota(label_w.begin(),label_w.end(),0);
	  label_w[order_w-1] = -1;
	  label_s[0] = -1;
	  std::iota(label_t.begin(),label_t.end(),0);
	  tcapi::contract(ctx,W[site_adrs],label_w,
			S,label_s,T_update[site_local_adrs],label_t);
	  V_update = tcapi::copy(ctx,T_update[site_local_adrs]);
          }
	}
#ifdef MCMC_DEBUG
	if( mpi_rank == 0 ) {
	  std::cout << " mcmc_update_line: end preparation of T_update " << std::endl;
	}
#endif
	// prepare
	IntT lenv_end_adrs = -1;
	IntT renv_end_adrs = static_cast<IntT>(line_target.size());
	if( sampling_direction == 0 ) {
	  lenv_end_adrs = local_adrs-1;
	} else {
	  renv_end_adrs = local_adrs+1;
	}
	std::vector<std::pair<IntT,IntT>> local_edges_temporary;
	std::vector<std::pair<IntT,IntT>> update_edges_temporary;
	mcmc_prep_local_target_tensors(ctx,line_target,edges,site_line_sets,
				       T_local,site_idx_local,
				       site_to_mpi_rank_local,
				       lenv_end_adrs,renv_end_adrs,
				       local_edges_temporary,comm);
	mcmc_prep_local_target_tensors(ctx,line_target,edges,site_line_sets,
				       T_update,site_idx_update,
				       site_to_mpi_rank_update,
				       lenv_end_adrs,renv_end_adrs,
				       update_edges_temporary,comm);

	if( sampling_direction == 0 ) {
	  if( local_adrs > min_local_adrs ) {
	    local_sites.push_back(site_label_nenv);
	    site_to_mpi_rank_local[site_label_nenv] = mpi_rank_nenv;
	    site_to_mpi_rank_update[site_label_nenv] = mpi_rank_nenv;
	    if( mpi_rank == mpi_rank_nenv ) {
	      T_local.push_back(tcapi::copy(ctx,Nenv));
	      site_idx_local.push_back(site_label_nenv);
	      T_update.push_back(tcapi::copy(ctx,Nenv));
	      site_idx_update.push_back(site_label_nenv);
	    }
	  }
	} else {
	  if( local_adrs < max_local_adrs ) {
	    local_sites.push_back(site_label_nenv);
	    site_to_mpi_rank_local[site_label_nenv] = mpi_rank_nenv;
	    site_to_mpi_rank_update[site_label_nenv] = mpi_rank_nenv;
	    if( mpi_rank == mpi_rank_nenv ) {
	      T_local.push_back(tcapi::copy(ctx,Nenv));
	      site_idx_local.push_back(site_label_nenv);
	      T_update.push_back(tcapi::copy(ctx,Nenv));
	      site_idx_update.push_back(site_label_nenv);
	    }
	  }
	}
	std::map<IntT,IntT> site_priority;
	std::map<IntT,IntT> label_map;
	site_priority[line_target[local_adrs]] = 1;
	auto edges_contract = extract_induced_edges(local_edges_temporary,
						    local_sites);
	auto sites_temporary = GetSiteIndexFromBond(local_edges_temporary);
	for(const auto & site : sites_temporary) {
	  auto it = site_to_mpi_rank_local.find(site);
	  if( it == site_to_mpi_rank_local.end() ) {
	    site_to_mpi_rank_local[site] = -1;
	    site_to_mpi_rank_update[site] = -1;
	  }
	}
#ifdef MCMC_DEBUG
	if( mpi_rank == 0 ) {
	  std::cout << " mcmc_update_line: start local contraction for sites ";
	  for(const auto & site : local_sites) {
	    std::cout << " " << site;
	  }
	  std::cout << ", edges =";
	  for(const auto & [u,v] : local_edges_temporary) {
	    std::cout << " (" << u << "," << v << ")";
	  }
	  std::cout << std::endl;
	}
#endif
	detail::graph_tensor_network_partial_contraction(ctx,local_edges_temporary,
					 T_local,site_idx_local,
					 site_to_mpi_rank_local,
					 edges_contract,
					 site_priority,
					 label_map,comm);
#ifdef MCMC_DEBUG
	if( mpi_rank == 0 ) {
	  std::cout << " mcmc_update_line: end local contraction, edges = ";
	  for(const auto & [u,v] : local_edges_temporary) {
	    std::cout << " (" << u << "," << v << ")";
	  }
	  std::cout << std::endl;
	  std::cout << " mcmc_update_line: start update contraction for edges =";
	  for(const auto & [u,v] : update_edges_temporary) {
	    std::cout << " (" << u << "," << v << ")";
	  }
	  std::cout << std::endl;
	}
#endif
	detail::graph_tensor_network_partial_contraction(ctx,update_edges_temporary,
					 T_update,site_idx_update,
					 site_to_mpi_rank_update,
					 edges_contract,
					 site_priority,
					 label_map,comm);
#ifdef MCMC_DEBUG
	if( mpi_rank == 0 ) {
	  std::cout << " mcmc_update_line: end update contraction, edges = ";
	  for(const auto & [u,v] : update_edges_temporary) {
	    std::cout << " (" << u << "," << v << ")";
	  }
	  std::cout << std::endl;
	}
#endif
	// keep new environment tensor
	TenT Nenv_update;
	site_label_nenv = line_target[local_adrs];
	mpi_rank_nenv = site_to_mpi_rank.at(site_label_nenv);
	if( mpi_rank == mpi_rank_nenv ) {
	  auto it_nenv_adrs = std::find(site_idx_local.begin(),
					site_idx_local.end(),
					site_label_nenv);
	  auto nenv_adrs = std::distance(site_idx_local.begin(),
					 it_nenv_adrs);
	  Nenv = tcapi::copy(ctx,T_local[nenv_adrs]);
	  Nenv_update = tcapi::copy(ctx,T_update[nenv_adrs]);
	}
#ifdef MCMC_DEBUG
	if( mpi_rank == 0 ) {
	  std::cout << " mcmc_update_line: keep target environment tensor as a candidate " << std::endl;
	}
#endif
	// contract tensor between new and remaining tensor
	TenT Result_local;
	TenT Result_update;
	if( local_point < static_cast<int>(site_line_sets.size())-1 ) {
	  // prepare the form of Nenv tensor
#ifdef MCMC_DEBUG
	  if( mpi_rank_nenv == mpi_rank ) {
	    auto it_adrs_nenv = std::find(site_idx_local.begin(),
					  site_idx_local.end(),
					  site_label_nenv);
	    auto adrs_nenv = std::distance(site_idx_local.begin(),
					   it_adrs_nenv);
	    auto shape_nenv = tcapi::shape(ctx,T_local[adrs_nenv]);
	    std::cout << " mcmc_update_line: start reshape of target local environment tensor = [";
	    for(const auto & dim : shape_nenv) {
	      std::cout << " " << dim;
	    }
	    std::cout << "]" << std::endl;
	  }
#endif
	  lenv_end_adrs = -1;
	  renv_end_adrs = static_cast<IntT>(line_target.size());
	  if( sampling_direction == 0 ) {
	    renv_end_adrs = local_adrs+1;
	  } else {
	    lenv_end_adrs = local_adrs-1;
	  }
	  auto edges_base = local_edges_temporary;
	  mcmc_prep_local_target_tensors(ctx,line_target,edges_base,
					 site_line_sets,
					 T_local,site_idx_local,
					 site_to_mpi_rank_local,
					 lenv_end_adrs,renv_end_adrs,
					 local_edges_temporary,comm);
	  mcmc_prep_local_target_tensors(ctx,line_target,edges_base,
					 site_line_sets,
					 T_update,site_idx_update,
					 site_to_mpi_rank_update,
					 lenv_end_adrs,renv_end_adrs,
					 update_edges_temporary,comm);
#ifdef MCMC_DEBUG
	  if( mpi_rank_nenv == mpi_rank ) {
	    auto it_adrs_nenv = std::find(site_idx_local.begin(),
					  site_idx_local.end(),
					  line_target[local_adrs]);
	    auto adrs_nenv = std::distance(site_idx_local.begin(),
					   it_adrs_nenv);
	    auto shape_nenv = tcapi::shape(ctx,T_local[adrs_nenv]);
	    std::cout << " mcmc_update_line: reshaped target local environment tensor = [";
	    for(const auto & dim : shape_nenv) {
	      std::cout << " " << dim;
	    }
	    std::cout << "]" << std::endl;
	  }
#endif
	  // prepare the form of Tenv tensor
	  edges_base = edges;
	  std::vector<TenT> T_reform;
	  std::vector<IntT> site_idx_reform;
	  std::vector<std::pair<IntT,IntT>> edges_reform;
	  std::map<IntT,IntT> site_to_mpi_rank_reform;
	  mcmc_prep_local_target_tensors(ctx,line_target,edges_base,
					 site_line_sets,
					 T_reform,site_idx_reform,
					 site_to_mpi_rank_reform,
					 lenv_end_adrs,renv_end_adrs,
					 edges_reform,comm);
	  IntT site_label_tenv = ( sampling_direction == 0 ) ?
	    line_target[local_adrs+1] : line_target[local_adrs-1];
	  auto mpi_rank_tenv = site_to_mpi_rank.at(site_label_tenv);
	  site_to_mpi_rank_reform[site_label_tenv] = mpi_rank_tenv;
	  if( mpi_rank_tenv == mpi_rank ) {
	    auto it_tenv_adrs = std::find(site_idx_tenv.begin(),
					  site_idx_tenv.end(),
					  site_label_tenv);
	    auto tenv_adrs = std::distance(site_idx_tenv.begin(),
					   it_tenv_adrs);
	    T_reform.push_back(tcapi::copy(ctx,T_tenv[tenv_adrs]));
	    site_idx_reform.push_back(site_label_tenv);
	  }
	  lenv_end_adrs = -1;
	  renv_end_adrs = static_cast<IntT>(line_target.size());
	  if( sampling_direction == 0 ) {
	    lenv_end_adrs = local_adrs;
	  } else {
	    renv_end_adrs = local_adrs;
	  }
	  edges_base = edges_reform;
#ifdef MCMC_DEBUG
	  if( mpi_rank_tenv == mpi_rank ) {
	    auto it_adrs_reform = std::find(site_idx_reform.begin(),
					  site_idx_reform.end(),
					  site_label_tenv);
	    auto adrs_reform = std::distance(site_idx_reform.begin(),
					     it_adrs_reform);
	    auto shape_reform = tcapi::shape(ctx,T_reform[adrs_reform]);
	    std::cout << " mcmc_update_line: start reshape of complementary environment tensor = [";
	    for(const auto & dim : shape_reform) {
	      std::cout << " " << dim;
	    }
	    std::cout << "]" << std::endl;
	  }
#endif
	  mcmc_prep_local_target_tensors(ctx,line_target,edges_base,
					 site_line_sets,
					 T_reform,site_idx_reform,
					 site_to_mpi_rank_reform,
					 lenv_end_adrs,renv_end_adrs,
					 edges_reform,comm);
#ifdef MCMC_DEBUG
	  if( mpi_rank_tenv == mpi_rank ) {
	    auto it_adrs_reform = std::find(site_idx_reform.begin(),
					  site_idx_reform.end(),
					  site_label_tenv);
	    auto adrs_reform = std::distance(site_idx_reform.begin(),
					     it_adrs_reform);
	    auto shape_reform = tcapi::shape(ctx,T_reform[adrs_reform]);
	    std::cout << " mcmc_update_line: reshaped complementary environment tensor = [";
	    for(const auto & dim : shape_reform) {
	      std::cout << " " << dim;
	    }
	    std::cout << "]" << std::endl;
	  }
#endif
	  // perform contraction between
	  int mpi_transfer_type = 0;
	  if( mpi_rank_tenv == mpi_rank ) {
	    mpi_transfer_type += 1;
	  }
	  if( mpi_rank_nenv == mpi_rank ) {
	    mpi_transfer_type += 2;
	  }
	  TenT Tenv;
	  if( mpi_transfer_type == 1 || mpi_transfer_type == 3 ) {
	    auto it_tenv_adrs = std::find(site_idx_reform.begin(),
					  site_idx_reform.end(),
					  site_label_tenv);
	    auto tenv_adrs = std::distance(site_idx_reform.begin(),
					   it_tenv_adrs);
	    Tenv = tcapi::copy(ctx,T_reform[tenv_adrs]);
	  }
	  if( mpi_transfer_type == 1 ) {
	    MpiSend(ctx,Tenv,mpi_rank_nenv,comm);
	  }
	  if( mpi_transfer_type == 2 ) {
	    MpiRecv(ctx,Tenv,mpi_rank_tenv,comm);
	  }
	  if( mpi_transfer_type == 2 || mpi_transfer_type == 3 ) {
	    auto it_nenv_adrs = std::find(site_idx_local.begin(),
					  site_idx_local.end(),
					  line_target[local_adrs]);
	    auto nenv_adrs = std::distance(site_idx_local.begin(),
					   it_nenv_adrs);
	    auto order_nenv = tcapi::order(ctx,T_local[nenv_adrs]);
	    auto shape_nenv = tcapi::shape(ctx,T_local[nenv_adrs]);
	    auto order_tenv = tcapi::order(ctx,Tenv);
	    auto shape_tenv = tcapi::shape(ctx,Tenv);
	    tcapi::List<BondLabelT> label_nenv(order_nenv);
	    tcapi::List<BondLabelT> label_tenv(order_tenv);
	    tcapi::List<BondLabelT> label_result;
	    if( order_nenv != order_tenv ) {
	      std::cerr << " mcmc_update_line: incompatible order of tensors between target and complementary environment tensors";
	      MPI_Abort(comm,1);
	    }
#ifdef MCMC_DEBUG
	    std::cout << " mcmc_update_line: start contraction of target and complementary environment tensors " << std::endl;
#endif
	    std::iota(label_nenv.begin(),label_nenv.end(),0);
	    std::iota(label_tenv.begin(),label_tenv.end(),0);
	    tcapi::contract(ctx,T_local[nenv_adrs],label_nenv,
			  Tenv,label_tenv,
			  Result_local,label_result);
	    tcapi::contract(ctx,T_update[nenv_adrs],label_nenv,
			  Tenv,label_tenv,
			  Result_update,label_result);
#ifdef MCMC_DEBUG
	    std::cout << " mcmc_update_line: end contraction of target and complementary environment tensors " << std::endl;
#endif
	  }
	} else {
	  // the case where no opposite side environment exists (because it is at the edge)
	  if( mpi_rank_nenv == mpi_rank ) {
	    auto it_site_local_adrs = std::find(site_idx_local.begin(),
						site_idx_local.end(),
						line_target[local_adrs]);
	    auto site_local_adrs = std::distance(site_idx_local.begin(),
						 it_site_local_adrs);
	    Result_local = tcapi::copy(ctx,T_local[site_local_adrs]);
	    Result_update = tcapi::copy(ctx,T_update[site_local_adrs]);
	  }
	}

        // Fixed sites keep Nenv from the current-state contraction and consume
        // neither proposal nor acceptance random numbers.
	if( mpi_rank_nenv == mpi_rank && !fixed_site ) {
	  auto order_r = tcapi::order(ctx,Result_local);
	  CoorsT coors_r(order_r,0);
	  ElemT res_local = tcapi::get_elem(ctx,Result_local,coors_r);
	  ElemT res_update = tcapi::get_elem(ctx,Result_update,coors_r);
	  RealT abs_local = std::abs(res_local);
	  RealT abs_update = std::abs(res_update);
	  RealT probability = sampling_method==0
            ? mcmc_pair_metropolis_probability(abs_local,abs_update)
            : mcmc_pair_heat_bath_probability(abs_local,abs_update);
	  RealT value = dist_real(rng);
	  if( value < probability ) {
	    Nenv = tcapi::copy(ctx,Nenv_update);
	    auto it_site_target_adrs = std::find(site_idx_target.begin(),
						 site_idx_target.end(),
						 line_target[local_adrs]);
	    auto site_target_adrs = std::distance(site_idx_target.begin(),
						  it_site_target_adrs);
	    T_target[site_target_adrs] = tcapi::copy(ctx,V_update);
	    auto it_site_adrs = std::find(site_idx.begin(),
					  site_idx.end(),
					  line_target[local_adrs]);
	    auto site_adrs = std::distance(site_idx.begin(),
					   it_site_adrs);
#ifdef MCMC_DEBUG
	    std::cout << " mcmc_update_line: update from state "
		      << site_states[site_adrs] << " to " << state_update
		      << " at site " << line_target[local_adrs]
		      << ", acceptance probability " << probability
		      << ", random number = " << value
		      << std::endl;
#endif
	    site_states[site_adrs] = state_update;
	  } else {
#ifdef MCMC_DEBUG
	    auto it_site_adrs = std::find(site_idx.begin(),
					  site_idx.end(),
					  line_target[local_adrs]);
	    auto site_adrs = std::distance(site_idx.begin(),
					   it_site_adrs);
	    std::cout << " mcmc_update_line: keep state "
		      << site_states[site_adrs]
		      << " at site " << line_target[local_adrs]
		      << ", acceptance probability " << probability
		      << ", random number = " << value
		      << std::endl;
#endif
	  }
	}
	// end update for each site update
      } // end for(int point_adrs (same with local_adrs))
    } // end single-site updates
#ifdef MCMC_DEBUG
    if( mpi_rank == 0 ) {
      std::cout << " mcmc_update_line: ================================================== end "
		<< std::endl;
    }
    MPI_Barrier(comm);
#endif
  }


  /**
   * @brief Initialize one chain using separately supplied output containers.
   * @tparam TenT TCAPI tensor type.
   * @tparam IntT Signed integer type for site labels, states and indices.
   * @param[in,out] ctx TCAPI tensor backend context.
   * @param[in] edges Ordered undirected graph edges; virtual axes follow incident-edge order.
   * @param[in] W Original rank-local tensors, with the physical axis last; never modified.
   * @param[in] site_idx Site labels in the same local order as the tensor/state arrays.
   * @param[in] site_to_mpi_rank Global site-to-owner-rank map; original ownership is unchanged.
   * @param[in] site_states Initial physical states in rank-local site_idx order.
   * @param[in] lines Ordered nonempty lines, from lower to upper boundary.
   * @param[out] V Attached tensors in local site_idx order.
   * @param[out] line_edges Reduced ordered graphs indexed by target line.
   * @param[out] T_upper Boundary tensors indexed [line][local tensor]; upper region absorbed onto that line.
   * @param[out] site_idx_upper Local site labels for each T_upper[line].
   * @param[in,out] E_upper BP messages indexed [line][local message], initialized where BP is run.
   * @param[in,out] edge_idx_upper Local BP edge indices per line, relative to that boundary induced graph.
   * @param[out] T_lower Boundary tensors indexed [line][local tensor]; lower region absorbed onto that line.
   * @param[out] site_idx_lower Local site labels for each T_lower[line].
   * @param[in,out] E_lower BP messages indexed [line][local message], initialized where BP is run.
   * @param[in,out] edge_idx_lower Local BP edge indices per line, relative to that boundary induced graph.
   * @param[out] site_line_sets Validated fixed grouping saved for subsequent sampling.
   * @param[in] max_bpsteps Maximum BP iterations per boundary truncation.
   * @param[in] bp_tolerance BP residual stopping tolerance.
   * @param[in] max_bond_dim Maximum retained boundary bond dimension.
   * @param[in] sv_min Singular-value threshold passed to boundary truncation.
   * @param[in] tg_err Target truncation error passed to boundary truncation.
   * @param[in] comm MPI communicator; all ranks must participate in the same call order.
   * @param[in] grouping Optional fixed three-line grouping. Null selects automatic grouping from the original graph.
   *
   * Choose groups from the original graph when grouping is null; otherwise validate
   * the supplied grouping. Initialize attachments, reduced line graphs and both
   * boundary caches, then save the same groups for all subsequent updates.
   * @pre W/state/owner arrays are consistent and lines form a valid sampling plan.
   * Collective on comm. Exceptions after initialization starts may leave partial outputs.
   * @see mcmc_init_attachment
   * @see mcmc_init_edges
   * @see mcmc_init_bmps
   */
  template <typename TenT, typename IntT>
  void mcmc_initialize(context_handle_t<TenT> & ctx,
       const std::vector<std::pair<IntT,IntT>> & edges,
       const std::vector<TenT> & W,
       const std::vector<IntT> & site_idx,
       const std::map<IntT,IntT> & site_to_mpi_rank,
       const std::vector<IntT> & site_states,
       const std::vector<std::vector<IntT>> & lines,
       std::vector<TenT> & V,
       std::vector<std::vector<std::pair<IntT,IntT>>> & line_edges,
       std::vector<std::vector<TenT>> & T_upper,
       std::vector<std::vector<IntT>> & site_idx_upper,
       std::vector<std::vector<TenT>> & E_upper,
       std::vector<std::vector<IntT>> & edge_idx_upper,
       std::vector<std::vector<TenT>> & T_lower,
       std::vector<std::vector<IntT>> & site_idx_lower,
       std::vector<std::vector<TenT>> & E_lower,
       std::vector<std::vector<IntT>> & edge_idx_lower,
       std::vector<std::vector<std::vector<IntT>>> & site_line_sets,
       IntT max_bpsteps, real_t<TenT> bp_tolerance,
       bond_dim_t<TenT> max_bond_dim, real_t<TenT> sv_min, real_t<TenT> tg_err,
       MPI_Comm comm,
       const MCMCGrouping<IntT>* grouping=nullptr) {

    MCMCGrouping<IntT> automatic;
    if(!grouping) { automatic=make_mcmc_grouping(edges,lines);grouping=&automatic; }
    else validate_mcmc_grouping(edges,lines,*grouping);
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);

    mcmc_init_attachment(ctx,edges,W,site_idx,site_to_mpi_rank,
			 site_states,lines,V);
#ifdef MCMC_DEBUG
    if( mpi_rank == 0 ) {
      std::cout << " mcmc_initialize: end mcmc_init_attachment at mcmc_initialize " << std::endl;
      for(std::size_t l=0; l < lines.size(); l++) {
	std::cout << " mcmc_initialize: line " << l << " =";
	for(const auto & site : lines[l]) {
	  std::cout << " " << site;
	}
	std::cout << std::endl;
      }
    }
#endif

    line_edges.resize(lines.size());
    std::size_t min_line_idx = 0;
    std::size_t max_line_idx = (lines.size() > 0) ? (lines.size()-1) : 0;
    for(std::size_t k=min_line_idx; k <= max_line_idx; k++) {
      mcmc_init_edges(edges,lines,static_cast<IntT>(k),line_edges[k],grouping);
#ifdef MCMC_DEBUG
      if( mpi_rank == 0 ) {
	std::cout << " mcmc_initialize: end mcmc_init_edges for line "
		  << k << " at mcmc_initialize " << std::endl;
	std::cout << " mcmc_initialize: edges for line "
		  << k << ":";
	for(const auto & [u,v] : line_edges[k]) {
	  std::cout << " (" << u << "," << v << ")";
	}
	std::cout << std::endl;
      }
#endif
    }
    mcmc_init_bmps(ctx,edges,
		   V,site_idx,site_to_mpi_rank,lines,
		   T_upper,site_idx_upper,
		   E_upper,edge_idx_upper,
		   T_lower,site_idx_lower,
		   E_lower,edge_idx_lower,
		   max_bpsteps,bp_tolerance,
		   max_bond_dim,sv_min,tg_err,
		   comm,grouping);
#ifdef MCMC_DEBUG
    if( mpi_rank == 0 ) {
      std::cout << " mcmc_initialize: end mcmc_init_bmps at mcmc_initialize " << std::endl;
    }
#endif

    site_line_sets=*grouping;

  }

  /**
   * @brief Generate samples by forward and reverse boundary-MPS sweeps.
   * @tparam TenT TCAPI tensor type.
   * @tparam IntT Signed integer type for site labels, states and indices.
   * @tparam RNG Uniform random bit generator type.
   * @param[in,out] ctx TCAPI tensor backend context.
   * @param[in] edges Ordered undirected graph edges; virtual axes follow incident-edge order.
   * @param[in] W Original rank-local tensors, with the physical axis last; never modified.
   * @param[in] site_idx Site labels in the same local order as the tensor/state arrays.
   * @param[in] site_to_mpi_rank Global site-to-owner-rank map; original ownership is unchanged.
   * @param[in] lines Ordered nonempty lines, from lower to upper boundary.
   * @param[in,out] sampling_states Preallocated [sample][batch][local site] states. Entry [0][batch] must match the initialized chain; each output sample is after a full round trip.
   * @param[in] line_edges Initialized reduced graph for each target line; unchanged.
   * @param[in,out] V Attached tensors indexed [batch][local site].
   * @param[in,out] T_upper Boundary tensors indexed [batch][line][local tensor].
   * @param[in] site_idx_upper Local boundary site labels indexed [line][local tensor], shared across batches.
   * @param[in,out] E_upper BP messages indexed [batch][line][local message].
   * @param[in,out] edge_idx_upper BP edge indices indexed [line][local message], refreshed during initialization of messages.
   * @param[in,out] T_lower Boundary tensors indexed [batch][line][local tensor].
   * @param[in] site_idx_lower Local boundary site labels indexed [line][local tensor], shared across batches.
   * @param[in,out] E_lower BP messages indexed [batch][line][local message].
   * @param[in,out] edge_idx_lower BP edge indices indexed [line][local message], refreshed during initialization of messages.
   * @param[in] site_line_sets Grouping saved at initialization; validated and reused unchanged.
   * @param[in] sampling_method 0 for Metropolis; 1 for Barker; 2 for full-state heat bath.
   * @param[in] update_sites Window size 1 or 2; size 2 requires method 2.
   * @param[in,out] rng Rank-local generator; method 2 draws on the leading target owner; methods 0/1 on site owners.
   * @param[in] max_bp_steps Maximum BP iterations per boundary truncation.
   * @param[in] bp_tolerance BP residual stopping tolerance.
   * @param[in] max_bond_dim Maximum retained boundary bond dimension.
   * @param[in] sv_min Singular-value threshold passed to boundary truncation.
   * @param[in] tg_err Target truncation error passed to boundary truncation.
   * @param[in] comm MPI communicator; all ranks must participate in the same call order.
   * @param[in] grouping Optional grouping used only to check agreement with site_line_sets; null still reuses the saved grouping.
   *
   * For each sample, sweep lines lower-to-upper with left-to-right site updates, then
   * upper-to-lower with right-to-left updates. Refresh the boundary behind the sweep
   * by BMPS absorption and BP/truncation, skipping the unused final boundary.
   * Copy the final state into the next sample slot to continue the chain. No burn-in
   * or thinning is performed internally. Empty sample storage performs no sweeps.
   * @pre All batches have compatible array sizes and initialized caches matching the
   * starting states. Use at least two nonempty lines; graph/plan metadata agree on ranks.
   * Collective on comm. Truncated boundary contraction can approximate the distribution.
   * @throws std::invalid_argument If the method or saved grouping is inconsistent.
   * @see mcmc_update_line
   */
  template <typename TenT, typename IntT, typename RNG>
  void mcmc_sampling(context_handle_t<TenT> & ctx,
       const std::vector<std::pair<IntT,IntT>> & edges,
       const std::vector<TenT> & W,
       const std::vector<IntT> & site_idx,
       const std::map<IntT,IntT> & site_to_mpi_rank,
       const std::vector<std::vector<IntT>> & lines,
       std::vector<std::vector<std::vector<IntT>>> & sampling_states,
       std::vector<std::vector<std::pair<IntT,IntT>>> & line_edges,
       std::vector<std::vector<TenT>> & V,
       std::vector<std::vector<std::vector<TenT>>> & T_upper,
       std::vector<std::vector<IntT>> & site_idx_upper,
       std::vector<std::vector<std::vector<TenT>>> & E_upper,
       std::vector<std::vector<IntT>> & edge_idx_upper,
       std::vector<std::vector<std::vector<TenT>>> & T_lower,
       std::vector<std::vector<IntT>> & site_idx_lower,
       std::vector<std::vector<std::vector<TenT>>> & E_lower,
       std::vector<std::vector<IntT>> & edge_idx_lower,
       std::vector<std::vector<std::vector<IntT>>> & site_line_sets,
       IntT sampling_method, RNG & rng,
       IntT max_bp_steps, real_t<TenT> bp_tolerance,
       bond_dim_t<TenT> max_bond_dim, real_t<TenT> sv_min, real_t<TenT> tg_err,
       MPI_Comm comm,
       const MCMCGrouping<IntT>* grouping=nullptr, int update_sites=1) {

    validate_mcmc_sampling_method(sampling_method);
    if(update_sites!=1 && update_sites!=2) throw std::invalid_argument("update_sites must be 1 or 2");
    if(update_sites==2 && sampling_method!=2) throw std::invalid_argument("two-site updates require method 2");
    // Reuse the initialized plan, including in automatic mode. Never regroup
    // from a contracted graph or independently for the opposite sweep direction.
    if(grouping && site_line_sets!=*grouping)
      throw std::invalid_argument("sampling groups differ from initialization");
    validate_mcmc_grouping(edges,lines,site_line_sets);
    grouping=&site_line_sets;
    using BondDimT = typename tcapi::tensor_traits<TenT>::bond_dim_t;
    using RealT = typename tcapi::tensor_traits<TenT>::real_t;

    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
#ifdef MCMC_DEBUG
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
#endif

    std::size_t num_samples = sampling_states.size();
    std::size_t num_batches = (sampling_states.size() > 0) ?
      sampling_states[0].size() : 0;
    std::size_t min_line_idx = 0;
    std::size_t max_line_idx = (lines.size()>0) ? (lines.size()-1) : 0;
    std::vector<TenT> T_dammy;
    std::vector<IntT> site_idx_dammy;
    std::vector<IntT> line_dammy;
    // Sampling-call lifetime: each owner set is created once, including both
    // sweep directions, and released after the final sample/batch.
    std::map<std::vector<int>,MPI_Comm> region_communicators;
    std::vector<std::array<mcmc_edge_line_plan<IntT>,2>> edge_plans(sampling_method==2?lines.size():0);
    try {
    if(sampling_method==2 && num_samples && num_batches) {
      // Initial input metadata only. Numerical boundary shapes are never shared
      // globally, and plan selection is independent of subsequent truncations.
      std::vector<unsigned long long> model(site_to_mpi_rank.size()+1,0);
      std::size_t k=0;
      for(const auto& entry:site_to_mpi_rank) {
        if(entry.second==mpi_rank) {
          const auto i=std::find(site_idx.begin(),site_idx.end(),entry.first)-site_idx.begin();
          const auto shape=tcapi::shape(ctx,W.at(i));
          if(shape.empty() || shape.back()<1 || shape.back()>std::numeric_limits<int>::max())
            throw std::invalid_argument("sampling model physical dimension");
          model[k]=shape.back();
          for(std::size_t axis=0;axis+1<shape.size();++axis)model.back()=std::max(model.back(),static_cast<unsigned long long>(shape[axis]));
        }
        ++k;
      }
      MPI_Allreduce(MPI_IN_PLACE,model.data(),static_cast<int>(model.size()),MPI_UNSIGNED_LONG_LONG,MPI_MAX,comm);
      std::map<IntT,int> counts;k=0;for(const auto& entry:site_to_mpi_rank)counts[entry.first]=static_cast<int>(model[k++]);
      std::set<std::vector<int>> owner_sets;
      for(std::size_t j=0;j<lines.size();++j)for(int direction=0;direction<2;++direction) {
        edge_plans[j][direction]=mcmc_build_edge_line_plan(line_edges[j],lines[j],site_line_sets[j],site_to_mpi_rank,
            direction,update_sites,model.back(),static_cast<std::uint64_t>(max_bond_dim),counts,
            mcmc_boundary_edge_model(line_edges[j],lines,j));
        const auto& sets=std::get<3>(edge_plans[j][direction]);owner_sets.insert(sets.begin(),sets.end());
      }
      mcmc_acquire_region_communicators(comm,0,
        std::vector<std::vector<int>>(owner_sets.begin(),owner_sets.end()),region_communicators);
    }
    // forward update
    for(std::size_t sample=0; sample < num_samples; sample++) {
      for(std::size_t batch=0; batch < num_batches; batch++) {
	std::vector<TenT> T_bmps;
	std::vector<IntT> site_idx_bmps;
	std::vector<std::pair<IntT,IntT>> edges_bmps(edges);
	std::map<IntT,IntT> site_to_mpi_rank_bmps(site_to_mpi_rank);
	for(std::size_t line_idx=min_line_idx;
	    line_idx <= max_line_idx; line_idx++) {
	  std::vector<TenT> T_target;
	  std::vector<IntT> site_idx_target;
	  for(std::size_t site_adrs=0;
	      site_adrs < lines[line_idx].size(); site_adrs++) {
	    auto site_i = lines[line_idx][site_adrs];
	    auto mpi_rank_i = site_to_mpi_rank.at(site_i);
	    if( mpi_rank == mpi_rank_i ) {
	      auto it_adrs_i = std::find(site_idx.begin(),
					 site_idx.end(),
					 site_i);
	      auto adrs_i = std::distance(site_idx.begin(),
					  it_adrs_i);
	      T_target.push_back(tcapi::copy(ctx,V[batch][adrs_i]));
	      site_idx_target.push_back(site_i);
	    }
	  }
	  if( line_idx == min_line_idx ) {
	    mcmc_update_line(ctx,line_edges[line_idx],
			     W,site_idx,site_to_mpi_rank,
			     lines[line_idx+1],T_upper[batch][line_idx+1],site_idx_upper[line_idx+1],
			     lines[line_idx],T_target,site_idx_target,
			     line_dammy,T_dammy,site_idx_dammy,
			     site_line_sets[line_idx],sampling_states[sample][batch],
			     comm,rng,sampling_method,0,update_sites,&region_communicators,sampling_method==2?&edge_plans[line_idx][0]:nullptr);
	  } else if ( line_idx == max_line_idx ) {
	    mcmc_update_line(ctx,line_edges[line_idx],
			     W,site_idx,site_to_mpi_rank,
			     line_dammy,T_dammy,site_idx_dammy,
			     lines[line_idx],T_target,site_idx_target,
			     lines[line_idx-1],T_lower[batch][line_idx-1],site_idx_lower[line_idx-1],
			     site_line_sets[line_idx],sampling_states[sample][batch],
			     comm,rng,sampling_method,0,update_sites,&region_communicators,sampling_method==2?&edge_plans[line_idx][0]:nullptr);
	  } else {
	    mcmc_update_line(ctx,line_edges[line_idx],
			     W,site_idx,site_to_mpi_rank,
			     lines[line_idx+1],T_upper[batch][line_idx+1],site_idx_upper[line_idx+1],
			     lines[line_idx],T_target,site_idx_target,
			     lines[line_idx-1],T_lower[batch][line_idx-1],site_idx_lower[line_idx-1],
			     site_line_sets[line_idx],sampling_states[sample][batch],
			     comm,rng,sampling_method,0,update_sites,&region_communicators,sampling_method==2?&edge_plans[line_idx][0]:nullptr);
	  }
	  // retrieve V and lower
	  for(std::size_t site_adrs=0;
	      site_adrs < lines[line_idx].size(); site_adrs++) {
	    auto site_i = lines[line_idx][site_adrs];
	    auto mpi_rank_i = site_to_mpi_rank.at(site_i);
	    if( mpi_rank == mpi_rank_i ) {
	      auto it_adrs_i = std::find(site_idx.begin(),
					 site_idx.end(),
					 site_i);
	      auto adrs_i = std::distance(site_idx.begin(),
					  it_adrs_i);
	      auto it_adrs_t = std::find(site_idx_target.begin(),
					 site_idx_target.end(),
					 site_i);
	      auto adrs_t = std::distance(site_idx_target.begin(),
					  it_adrs_t);
	      V[batch][adrs_i] = tcapi::copy(ctx,T_target[adrs_t]);
	    }
	  }
	  if( line_idx < max_line_idx ) {
	    for(std::size_t site_adrs=0;
		site_adrs < lines[line_idx].size(); site_adrs++) {
	      auto site_i = lines[line_idx][site_adrs];
	      auto mpi_rank_i = site_to_mpi_rank.at(site_i);
	      if( mpi_rank == mpi_rank_i ) {
		auto it_adrs_i = std::find(site_idx.begin(),
					   site_idx.end(),
					   site_i);
		auto adrs_i = std::distance(site_idx.begin(),
					    it_adrs_i);
		T_bmps.push_back(tcapi::copy(ctx,V[batch][adrs_i]));
		site_idx_bmps.push_back(site_i);
	      }
	    }
	    if( line_idx > min_line_idx ) {
	      // perform bmps incrementation
#ifdef MCMC_DEBUG
	      if( mpi_rank == 0 ) {
		std::cout << " mcmc_sampling: start boundary-mps increment of lines "
			  << line_idx-1 << " and " << line_idx << std::endl;
		std::cout << " mcmc_sampling: base site information =";
		for(const auto & site : lines[line_idx-1]) {
		  std::cout << " " << site << " (rank " << site_to_mpi_rank_bmps.at(site) << ")";
		}
		std::cout << std::endl;
		std::cout << " mcmc_sampling: attached site information =";
		for(const auto & site : lines[line_idx]) {
		  std::cout << " " << site << " (rank " << site_to_mpi_rank_bmps.at(site) << ")";
		}
		std::cout << std::endl;
	      }
#endif
	      bmps_increment_line(ctx,
				  lines[line_idx-1],lines[line_idx],
				  edges_bmps,
				  T_bmps,site_idx_bmps,
				  site_to_mpi_rank_bmps,
				  comm, &(*grouping)[line_idx]);
	      // perform belief propagation and truncation
	      std::vector<BondDimT> res_bond_dim;
	      std::vector<RealT> res_trunc_err;
#ifdef MCMC_DEBUG
	      if( mpi_rank == 0 ) {
		std::cout << " mcmc_sampling: start boundary-mps truncation of line "
			  << line_idx << std::endl;
		std::cout << " mcmc_sampling: site information =";
		for(const auto & site : lines[line_idx]) {
		  std::cout << " " << site << " (rank " << site_to_mpi_rank_bmps.at(site) << ")";
		}
		std::cout << std::endl;
	      }
#endif
	      bmps_bp_truncation(ctx,edges_bmps,lines[line_idx],
				 T_bmps,site_idx_bmps,site_to_mpi_rank_bmps,
				 E_lower[batch][line_idx],edge_idx_lower[line_idx],
				 max_bp_steps,bp_tolerance,
				 max_bond_dim,sv_min,tg_err,res_bond_dim,
				 res_trunc_err,comm,true);
#ifdef MCMC_DEBUG
	      if( mpi_rank == 0 ) {
		std::cout << " mcmc_sampling: end boundary-mps trunctaion of line "
			  << line_idx << std::endl;
	      }
#endif
	    }
	    // update T_lower[batch][line_idx]
	    for(std::size_t site_adrs=0;
		site_adrs < lines[line_idx].size(); site_adrs++) {
	      auto site_i = lines[line_idx][site_adrs];
	      auto mpi_rank_i = site_to_mpi_rank.at(site_i);
	      if( mpi_rank == mpi_rank_i ) {
		auto it_adrs_bmps = std::find(site_idx_bmps.begin(),
					      site_idx_bmps.end(),
					      site_i);
		auto adrs_bmps = std::distance(site_idx_bmps.begin(),
					       it_adrs_bmps);
		auto it_adrs_lower = std::find(site_idx_lower[line_idx].begin(),
					       site_idx_lower[line_idx].end(),
					       site_i);
		auto adrs_lower = std::distance(site_idx_lower[line_idx].begin(),
						it_adrs_lower);
		T_lower[batch][line_idx][adrs_lower] = tcapi::copy(ctx,T_bmps[adrs_bmps]);
	      }
	    }
	  } // end if(line_idx < max_line_idx)
	} // end for(line_idx)
      } // end for(batch)

#ifdef MCMC_DEBUG
      if( mpi_rank == 0 ) {
	std::cout << " mcmc_sampling: ---- end lower to upper updates ----- " << std::endl;
	std::cout << " mcmc_sampling: start upper to lower updates ----- " << std::endl;
      }
#endif

      // from upper to lower
      for(std::size_t batch=0; batch < num_batches; batch++) {
	std::vector<TenT> T_bmps;
	std::vector<IntT> site_idx_bmps;
	std::vector<std::pair<IntT,IntT>> edges_bmps(edges);
	std::map<IntT,IntT> site_to_mpi_rank_bmps(site_to_mpi_rank);
	for(std::size_t rline_idx=0;
	    rline_idx < lines.size(); rline_idx++) {
	  std::size_t line_idx = max_line_idx - rline_idx;
	  std::vector<TenT> T_target;
	  std::vector<IntT> site_idx_target;
	  for(std::size_t site_adrs=0;
	      site_adrs < lines[line_idx].size(); site_adrs++) {
	    auto site_i = lines[line_idx][site_adrs];
	    auto mpi_rank_i = site_to_mpi_rank.at(site_i);
	    if( mpi_rank == mpi_rank_i ) {
	      auto it_adrs_i = std::find(site_idx.begin(),
					 site_idx.end(),
					 site_i);
	      auto adrs_i = std::distance(site_idx.begin(),it_adrs_i);
	      T_target.push_back(tcapi::copy(ctx,V[batch][adrs_i]));
	      site_idx_target.push_back(site_i);
	    }
	  } // end for(site_adrs)
	  if( line_idx == min_line_idx ) {
	    mcmc_update_line(ctx,line_edges[line_idx],
			     W,site_idx,site_to_mpi_rank,
			     lines[line_idx+1],T_upper[batch][line_idx+1],site_idx_upper[line_idx+1],
			     lines[line_idx],T_target,site_idx_target,
			     line_dammy,T_dammy,site_idx_dammy,
			     site_line_sets[line_idx],sampling_states[sample][batch],
			     comm,rng,sampling_method,1,update_sites,&region_communicators,sampling_method==2?&edge_plans[line_idx][1]:nullptr);
	  } else if ( line_idx == max_line_idx ) {
	    mcmc_update_line(ctx,line_edges[line_idx],
			     W,site_idx,site_to_mpi_rank,
			     line_dammy,T_dammy,site_idx_dammy,
			     lines[line_idx],T_target,site_idx_target,
			     lines[line_idx-1],T_lower[batch][line_idx-1],site_idx_lower[line_idx-1],
			     site_line_sets[line_idx],sampling_states[sample][batch],
			     comm,rng,sampling_method,1,update_sites,&region_communicators,sampling_method==2?&edge_plans[line_idx][1]:nullptr);
	  } else {
	    mcmc_update_line(ctx,line_edges[line_idx],
			     W,site_idx,site_to_mpi_rank,
			     lines[line_idx+1],T_upper[batch][line_idx+1],site_idx_upper[line_idx+1],
			     lines[line_idx],T_target,site_idx_target,
			     lines[line_idx-1],T_lower[batch][line_idx-1],site_idx_lower[line_idx-1],
			     site_line_sets[line_idx],sampling_states[sample][batch],
			     comm,rng,sampling_method,1,update_sites,&region_communicators,sampling_method==2?&edge_plans[line_idx][1]:nullptr);
	  }
	  // retrieve V
	  for(std::size_t site_adrs=0; site_adrs < lines[line_idx].size(); site_adrs++) {
	    auto site_i = lines[line_idx][site_adrs];
	    auto mpi_rank_i = site_to_mpi_rank.at(site_i);
	    if( mpi_rank == mpi_rank_i ) {
	      auto it_adrs_i = std::find(site_idx.begin(),site_idx.end(),
					 site_i);
	      auto adrs_i = std::distance(site_idx.begin(),it_adrs_i);
	      auto it_adrs_t = std::find(site_idx_target.begin(),site_idx_target.end(),
					 site_i);
	      auto adrs_t = std::distance(site_idx_target.begin(),
					  it_adrs_t);
	      V[batch][adrs_i] = tcapi::copy(ctx,T_target[adrs_t]);
	    }
	  }
	  if( line_idx > min_line_idx ) {
	    for(std::size_t site_adrs=0; site_adrs < lines[line_idx].size(); site_adrs++) {
	      auto site_i = lines[line_idx][site_adrs];
	      auto mpi_rank_i = site_to_mpi_rank.at(site_i);
	      if( mpi_rank == mpi_rank_i ) {
		auto it_adrs_i = std::find(site_idx.begin(),
					   site_idx.end(),
					   site_i);
		auto adrs_i = std::distance(site_idx.begin(),it_adrs_i);
		T_bmps.push_back(tcapi::copy(ctx,V[batch][adrs_i]));
		site_idx_bmps.push_back(site_i);
	      }
	    }
	    if( line_idx < max_line_idx ) {
	      // perform bmps incrementation
#ifdef MCMC_DEBUG
	      if( mpi_rank == 0 ) {
		std::cout << " mcmc_sampling: start boundary-mps increment of lines "
			  << line_idx+1 << " and " << line_idx << std::endl;
		std::cout << " mcmc_sampling: base site information =";
		for(const auto & site : lines[line_idx+1]) {
		  std::cout << " " << site << " (rank " << site_to_mpi_rank_bmps.at(site) << ")";
		}
		std::cout << std::endl;
		std::cout << " mcmc_sampling: attached site information =";
		for(const auto & site : lines[line_idx]) {
		  std::cout << " " << site << " (rank " << site_to_mpi_rank_bmps.at(site) << ")";
		}
		std::cout << std::endl;
	      }
#endif
	      bmps_increment_line(ctx,lines[line_idx+1],lines[line_idx],
				  edges_bmps,T_bmps,site_idx_bmps,
				  site_to_mpi_rank_bmps,
				  comm, &(*grouping)[line_idx]);
	      // perform belief propagation and truncation
	      std::vector<BondDimT> res_bond_dim;
	      std::vector<RealT> res_trunc_err;
#ifdef MCMC_DEBUG
	      if( mpi_rank == 0 ) {
		std::cout << " mcmc_sampling: start boundary-mps truncation of line "
			  << line_idx << std::endl;
		std::cout << " mcmc_sampling: site information =";
		for(const auto & site : lines[line_idx]) {
		  std::cout << " " << site << " (rank " << site_to_mpi_rank_bmps.at(site) << ")";
		}
		std::cout << std::endl;
	      }
#endif
	      bmps_bp_truncation(ctx,edges_bmps,lines[line_idx],
				 T_bmps,site_idx_bmps,site_to_mpi_rank_bmps,
				 E_upper[batch][line_idx],edge_idx_upper[line_idx],
				 max_bp_steps,bp_tolerance,
				 max_bond_dim,sv_min,tg_err,
				 res_bond_dim,res_trunc_err,
				 comm,true);
#ifdef MCMC_DEBUG
	      if( mpi_rank == 0 ) {
		std::cout << " mcmc_sampling: end boundary-mps trunctaion of line "
			  << line_idx << std::endl;
	      }
#endif
	    }
	    // update T_upper[batch][line_idx]
	    for(std::size_t site_adrs=0;
		site_adrs < lines[line_idx].size(); site_adrs++) {
	      auto site_i = lines[line_idx][site_adrs];
	      auto mpi_rank_i = site_to_mpi_rank.at(site_i);
	      if( mpi_rank == mpi_rank_i ) {
		auto it_adrs_bmps = std::find(site_idx_bmps.begin(),
					      site_idx_bmps.end(),
					      site_i);
		auto adrs_bmps = std::distance(site_idx_bmps.begin(),
					       it_adrs_bmps);
		auto it_adrs_upper = std::find(site_idx_upper[line_idx].begin(),
					       site_idx_upper[line_idx].end(),
					       site_i);
		auto adrs_upper = std::distance(site_idx_upper[line_idx].begin(),
						it_adrs_upper);
		T_upper[batch][line_idx][adrs_upper] = tcapi::copy(ctx,T_bmps[adrs_bmps]);
	      }
	    }
	  }
	} // end for(line_idx)

	if( sample < num_samples-1 ) {
	  sampling_states[sample+1][batch] = sampling_states[sample][batch];
	}

      } // end for(batch)
    } // end for(sample)
    mcmc_release_region_communicators(region_communicators);
    } catch(const std::exception& e) {
      if(sampling_method==2) detail::abort_tensor_operation(comm,e.what());
      throw;
    } catch(...) {
      if(sampling_method==2) detail::abort_tensor_operation(comm,"unknown sampling failure");
      throw;
    }
  }

  /**
   * @brief Persistent batched sampling caches and the shared contraction plan.
   * @tparam TenT TCAPI tensor type.
   * @tparam IntT Signed integer type for site labels, states and indices.
   *
   * Set lines and optionally manual_grouping before mcmc_initialize(). Tensor
   * storage is local to each MPI rank; topology and grouping agree across ranks.
   * Reinitialize after changing the graph, original tensors or grouping.
   * BP messages E_upper/E_lower are distinct from temporary left/right environments.
   * @see mcmc_initialize
   * @see mcmc_sampling
   */
  template <typename TenT, typename IntT>
  struct MCMCWorkspace {
    std::vector<std::vector<IntT>> lines; ///< Ordered site lists [line][position], lower to upper.
    MCMCGrouping<IntT> manual_grouping; ///< Optional input grouping; empty selects automatic planning.
    std::vector<std::vector<std::pair<IntT,IntT>>> line_edges; ///< Shared reduced graphs [target line][edge].
    std::vector<std::vector<TenT>> V; ///< Attached local tensors [batch][local site].
    std::vector<std::vector<std::vector<TenT>>> T_upper; ///< Upper boundary tensors [batch][line][local tensor].
    std::vector<std::vector<IntT>> site_idx_upper; ///< Boundary labels [line][local tensor], shared across batches.
    std::vector<std::vector<std::vector<TenT>>> E_upper; ///< Upper boundary BP messages [batch][line][local message].
    std::vector<std::vector<IntT>> edge_idx_upper; ///< BP edge indices [line][local message], shared across batches.
    std::vector<std::vector<std::vector<TenT>>> T_lower; ///< Lower boundary tensors [batch][line][local tensor].
    std::vector<std::vector<IntT>> site_idx_lower; ///< Boundary labels [line][local tensor], shared across batches.
    std::vector<std::vector<std::vector<TenT>>> E_lower; ///< Lower boundary BP messages [batch][line][local message].
    std::vector<std::vector<IntT>> edge_idx_lower; ///< BP edge indices [line][local message], shared across batches.
    std::vector<std::vector<std::vector<IntT>>> site_line_sets; ///< Resolved grouping [line][target position][group member], reused in both directions.
  };

  /**
   * @brief Initialize batched chains in an MCMCWorkspace.
   * @tparam TenT TCAPI tensor type.
   * @tparam IntT Signed integer type for site labels, states and indices.
   * @param[in,out] ctx TCAPI tensor backend context.
   * @param[in] edges Ordered undirected graph edges; virtual axes follow incident-edge order.
   * @param[in] W Original rank-local tensors, with the physical axis last; never modified.
   * @param[in] site_idx Site labels in the same local order as the tensor/state arrays.
   * @param[in] site_to_mpi_rank Global site-to-owner-rank map; original ownership is unchanged.
   * @param[in] site_states Initial states indexed [batch][local site], matching W/site_idx.
   * @param[in] max_bpsteps Maximum BP iterations per boundary truncation.
   * @param[in] bp_tolerance BP residual stopping tolerance.
   * @param[in] max_bond_dim Maximum retained boundary bond dimension.
   * @param[in] sv_min Singular-value threshold passed to boundary truncation.
   * @param[in] tg_err Target truncation error passed to boundary truncation.
   * @param[in,out] work Set lines and optional manual_grouping before calling; receives per-batch caches and shared graph metadata.
   * @param[in] comm MPI communicator; all ranks must participate in the same call order.
   *
   * Resolve a single grouping before processing batches. An empty manual_grouping
   * selects automatic planning; the resolved plan is stored in site_line_sets.
   * Each batch calls the individual-container initialization overload. Shared labels
   * and graphs are taken from the first batch. With zero batches, groups are validated
   * but per-line metadata is not rebuilt; this is not a populated sampling workspace.
   * Collective on comm; local states must be valid for the physical dimensions.
   * @see MCMCWorkspace
   */
  template <typename TenT, typename IntT>
  void mcmc_initialize(context_handle_t<TenT> & ctx,
		       const std::vector<std::pair<IntT,IntT>> & edges,
		       const std::vector<TenT> & W,
		       const std::vector<IntT> & site_idx,
		       const std::map<IntT,IntT> & site_to_mpi_rank,
		       const std::vector<std::vector<IntT>> & site_states,
		       IntT max_bpsteps, real_t<TenT> bp_tolerance,
		       bond_dim_t<TenT> max_bond_dim,
		       real_t<TenT> sv_min,
		       real_t<TenT> tg_err,
		       MCMCWorkspace<TenT,IntT> & work,
		       MPI_Comm comm) {

    // Validate even an empty batch request, before changing workspace outputs.
    const auto plan=work.manual_grouping.empty() ? make_mcmc_grouping(edges,work.lines) : work.manual_grouping;
    if(!work.manual_grouping.empty()) validate_mcmc_grouping(edges,work.lines,plan);
    work.site_line_sets=plan;
    work.V.resize(site_states.size());
    work.T_upper.resize(site_states.size());
    work.E_upper.resize(site_states.size());
    work.T_lower.resize(site_states.size());
    work.E_lower.resize(site_states.size());
    for(std::size_t batch=0; batch < site_states.size(); batch++) {
      auto lines=work.lines;
      std::vector<std::vector<std::pair<IntT,IntT>>> line_edges;
      std::vector<std::vector<IntT>> site_idx_upper;
      std::vector<std::vector<IntT>> edge_idx_upper;
      std::vector<std::vector<IntT>> site_idx_lower;
      std::vector<std::vector<IntT>> edge_idx_lower;
      std::vector<std::vector<std::vector<IntT>>> site_line_sets;
      mcmc_initialize(ctx,edges,W,site_idx,site_to_mpi_rank,
		      site_states[batch],
		      lines,work.V[batch],line_edges,
		      work.T_upper[batch],site_idx_upper,
		      work.E_upper[batch],edge_idx_upper,
		      work.T_lower[batch],site_idx_lower,
		      work.E_lower[batch],edge_idx_lower,
		      site_line_sets,
		      max_bpsteps,bp_tolerance,max_bond_dim,
		      sv_min,tg_err,comm,&plan);
      if( batch == 0 ) {
	work.lines = lines;
	work.line_edges = line_edges;
	work.site_idx_upper = site_idx_upper;
	work.edge_idx_upper = edge_idx_upper;
	work.site_idx_lower = site_idx_lower;
	work.edge_idx_lower = edge_idx_lower;
	work.site_line_sets = site_line_sets;
      }
    }
  }

  /**
   * @brief Sample batched chains using an initialized MCMCWorkspace.
   * @tparam TenT TCAPI tensor type.
   * @tparam IntT Signed integer type for site labels, states and indices.
   * @tparam RNG Uniform random bit generator type.
   * @param[in,out] ctx TCAPI tensor backend context.
   * @param[in] edges Ordered undirected graph edges; virtual axes follow incident-edge order.
   * @param[in] W Original rank-local tensors, with the physical axis last; never modified.
   * @param[in] site_idx Site labels in the same local order as the tensor/state arrays.
   * @param[in] site_to_mpi_rank Global site-to-owner-rank map; original ownership is unchanged.
   * @param[in,out] sampling_states Preallocated [sample][batch][local site] states. Entry [0][batch] must match the initialized chain; each output sample is after a full round trip.
   * @param[in] sampling_method 0 for Metropolis; 1 for Barker; 2 for full-state heat bath.
   * @param[in] update_sites Window size 1 or 2; size 2 requires method 2.
   * @param[in,out] rng Rank-local generator; method 2 draws on the leading target owner; methods 0/1 on site owners.
   * @param[in] max_bp_steps Maximum BP iterations per boundary truncation.
   * @param[in] bp_tolerance BP residual stopping tolerance.
   * @param[in] max_bond_dim Maximum retained boundary bond dimension.
   * @param[in] sv_min Singular-value threshold passed to boundary truncation.
   * @param[in] tg_err Target truncation error passed to boundary truncation.
   * @param[in,out] work Initialized workspace for the same graph, tensors, lines and starting states.
   * @param[in] comm MPI communicator; all ranks must participate in the same call order.
   *
   * Forwards workspace members to the individual-container sampling overload.
   * The saved site_line_sets is reused even in automatic mode. Changing the manual
   * grouping after initialization is rejected if it differs from that saved plan.
   * For each sample, sweep lines lower-to-upper with left-to-right site updates, then
   * upper-to-lower with right-to-left updates. Refresh the boundary behind the sweep
   * by BMPS absorption and BP/truncation, skipping the unused final boundary.
   * Copy the final state into the next sample slot to continue the chain. No burn-in
   * or thinning is performed internally. Empty sample storage performs no sweeps.
   * @pre All batches have compatible array sizes and initialized caches matching the
   * starting states. Use at least two nonempty lines; graph/plan metadata agree on ranks.
   * Collective on comm. Truncated boundary contraction can approximate the distribution.
   * @throws std::invalid_argument If the method or saved grouping is inconsistent.
   * @see mcmc_update_line
   */
  template <typename TenT, typename IntT, typename RNG>
  void mcmc_sampling(context_handle_t<TenT> & ctx,
       const std::vector<std::pair<IntT,IntT>> & edges,
       const std::vector<TenT> & W,
       const std::vector<IntT> & site_idx,
       const std::map<IntT,IntT> & site_to_mpi_rank,
       std::vector<std::vector<std::vector<IntT>>> & sampling_states,
       IntT sampling_method, RNG & rng,
       IntT max_bp_steps, real_t<TenT> bp_tolerance,
       bond_dim_t<TenT> max_bond_dim, real_t<TenT> sv_min, real_t<TenT> tg_err,
       MCMCWorkspace<TenT,IntT> & work,
       MPI_Comm comm, int update_sites=1) {
    mcmc_sampling(ctx,edges,W,site_idx,site_to_mpi_rank,work.lines,
		  sampling_states,work.line_edges,work.V,
		  work.T_upper,work.site_idx_upper,
		  work.E_upper,work.edge_idx_upper,
		  work.T_lower,work.site_idx_lower,
		  work.E_lower,work.edge_idx_lower,
		  work.site_line_sets,
		  sampling_method,rng,max_bp_steps,bp_tolerance,
		  max_bond_dim,sv_min,tg_err,comm,
                  work.manual_grouping.empty() ? nullptr : &work.manual_grouping,update_sites);
  }


  /** @brief Sample with full-state roulette heat bath by default.
   * Window size is 1 unless update_sites=2; all local states are enumerated.
   * Other methods remain available through the explicit-method overload.
   */
  template <typename TenT, typename IntT, typename RNG>
  void mcmc_sampling(context_handle_t<TenT>& ctx,
       const std::vector<std::pair<IntT,IntT>>& edges,
       const std::vector<TenT>& W, const std::vector<IntT>& site_idx,
       const std::map<IntT,IntT>& site_to_mpi_rank,
       std::vector<std::vector<std::vector<IntT>>>& sampling_states,
       RNG& rng, IntT max_bp_steps, real_t<TenT> bp_tolerance,
       bond_dim_t<TenT> max_bond_dim, real_t<TenT> sv_min, real_t<TenT> tg_err,
       MCMCWorkspace<TenT,IntT>& work, MPI_Comm comm, int update_sites=1) {
    mcmc_sampling(ctx,edges,W,site_idx,site_to_mpi_rank,sampling_states,IntT(2),rng,
                  max_bp_steps,bp_tolerance,max_bond_dim,sv_min,tg_err,work,comm,update_sites);
  }



}

#endif
