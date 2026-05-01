/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/ptns/function/bmps.h
@brief Functions for performing boundary mps
*/
#ifndef TNBP_PTNS_FUNCTION_BMPS_H
#define TNBP_PTNS_FUNCTION_BMPS_H

#include "tnbp/framework/typedef.h"
#include "tnbp/framework/graph.h"
#include "tnbp/framework/helper.h"
#include "tnbp/framework/mpiutility.h"

namespace tnbp {

  template <typename IntT>
  void bmps_setup_lines(const std::vector<std::pair<IntT,IntT>> & edges,
			const std::vector<IntT> & initial_boundary,
			std::vector<std::vector<IntT>> & lines) {
    lines = BuildContractionLayers(edges,initial_boundary);
  }

  template <typename TenT, typename IntT>
  void bmps_increment_line(context_handle_t<TenT> & ctx,
			   const std::vector<IntT> & line_in,
			   const std::vector<IntT> & line_out,
			   std::vector<std::pair<IntT,IntT>> & edges,
			   std::vector<TenT> & T,
			   std::vector<IntT> & site_idx,
			   std::map<IntT,IntT> & site_to_mpi_rank,
			   MPI_Comm comm) {
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    std::size_t min_adrs_line_in = 0;
    std::size_t max_adrs_line_in = (line_in.size()>0) ? (line_in.size()-1) : 0;
    for(std::size_t adrs_line_in = min_adrs_line_in;
	adrs_line_in <= max_adrs_line_in;
	++adrs_line_in) {
      auto site_in = line_in[adrs_line_in];
      auto bond_in = GetSurroundingBondIndex(site_in,edges);
      std::vector<IntT> bond_type(bond_in.size(),2);
      for(size_t k=0; k < bond_in.size(); k++) {
	auto site_a = edges[bond_in[k]].first;
	auto site_b = edges[bond_in[k]].second;
	auto site_next = (site_a == site_in) ? site_b : site_a;
	auto it_is_zero = std::find(line_out.begin(),
				    line_out.end(),
				    site_next);
	auto it_is_one  = std::find(line_in.begin(),
				    line_in.end(),
				    site_next);
	if( it_is_zero != line_out.end() ) {
	  bond_type[k] = 0;
	}
	if( it_is_one  != line_in.end() ) {
	  bond_type[k] = 1;
	}
      }
      auto rit_type_zero = std::find(bond_type.rbegin(),
				    bond_type.rend(),
				    0);
      auto it_type_one  = std::find(bond_type.begin(),
				    bond_type.end(),
				    1);
      size_t target_bond_index;
      IntT target_edge_adrs;
      if( rit_type_zero != bond_type.rend() ) {
	auto it_type_zero = rit_type_zero.base()-1;
	auto target_bond_index = std::distance(bond_type.begin(),
					       it_type_zero);
	target_edge_adrs = bond_in[target_bond_index];
      } else if (it_type_one != bond_type.end()) {
	auto target_bond_index = std::distance(bond_type.begin(),
					       it_type_one);
	target_edge_adrs = bond_in[target_bond_index];
      }
      auto site_a = edges[target_edge_adrs].first;
      auto site_b = edges[target_edge_adrs].second;
      auto site_next = (site_a == site_in) ? site_b : site_a;
      graph_tensor_contraction(ctx,edges,T,site_idx,site_to_mpi_rank,
			       site_next,site_in,site_next,comm);
#ifdef BMPS_DEBUG
      std::cout << " edges in bmps_increment_line:";
      for(const auto & [u,v] : edges) {
	std::cout << " (" << u << "," << v << ")";
      }
      std::cout << std::endl;
#endif
    }
  }
  
  template <typename TenT, typename IntT>
  void bmps_increment_line_for_edges(const std::vector<IntT> & line_in,
				     const std::vector<IntT> & line_out,
				     std::vector<std::pair<IntT,IntT>> & edges) {
    for(size_t adrs_line_in=0;
	adrs_line_in < line_in.size();
	++adrs_line_in) {
      auto site_in = line_in[adrs_line_in];
      auto bond_in = GetSurroundingBondIndex(site_in,edges);
      std::vector<IntT> bond_type(bond_in.size(),2);
      for(size_t k=0; k < bond_in.size(); k++) {
	auto site_a = edges[bond_in[k]].first;
	auto site_b = edges[bond_in[k]].second;
	auto site_next = (site_a == site_in) ? site_b : site_a;
	auto it_is_zero = std::find(line_out.begin(),
				    line_out.end(),
				    site_next);
	auto it_is_one  = std::find(line_in.begin(),
				    line_in.end(),
				    site_next);
	if( it_is_zero != line_out.end() ) {
	  bond_type[k] = 0;
	}
	if( it_is_one  != line_in.end() ) {
	  bond_type[k] = 1;
	}
      }
      auto it_type_zero = std::find(bond_type.begin(),
				    bond_type.end(),
				    0);
      auto it_type_one  = std::find(bond_type.begin(),
				    bond_type.end(),
				    1);
      size_t target_bond_index;
      IntT target_edge_adrs;
      if( it_type_zero != bond_type.end() ) {
	auto target_bond_index = std::distance(bond_type.begin(),
					       it_type_zero);
	target_edge_adrs = bond_in[target_bond_index];
      } else if (it_type_one != bond_type.end()) {
	auto target_bond_index = std::distance(bond_type.begin(),
					       it_type_one);
	target_edge_adrs = bond_in[target_bond_index];
      }
      auto site_a = edges[target_edge_adrs].first;
      auto site_b = edges[target_edge_adrs].second;
      auto site_next = (site_a == site_in) ? site_b : site_a;
      graph_tensor_contraction_for_edges(edges,
			       site_next,site_in,site_next);
#ifdef BMPS_DEBUG
      std::cout << " edges in bmps_increment_line_for_edges:";
      for(const auto & [u,v] : edges) {
	std::cout << " (" << u << "," << v << ")";
      }
      std::cout << std::endl;
#endif
    }
  }

  // Function to perform bp and truncation
  template <typename TenT, typename IntT, typename StepT>
  void bmps_bp_truncation(context_handle_t<TenT> ctx,
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
			  MPI_Comm comm,
			  bool do_init = true) {
    
    using RealT = typename tci::tensor_traits<TenT>::real_t;
    using BondDimT = typename tci::tensor_traits<TenT>::bond_dim_t;
    
    auto edges_for_bp = extract_induced_edges(edges,lines);
    if( do_init ) {
      init_edge_messenger_tensors(ctx,edges_for_bp,T,
				  site_idx,site_to_mpi_rank,
				  E,edge_idx,comm);
    }
    auto layering_edges = GreedyEdgeLayering(edges_for_bp);
    for(StepT step=0; step < max_bp_steps; step++) {
      std::vector<TenT> E_new;
      for(std::size_t layer=0; layer < layering_edges.size(); layer++) {
	BeliefPropagation(ctx,edges_for_bp,
			  T,site_idx,site_to_mpi_rank,
			  layering_edges[layer],
			  E,edge_idx,comm,E_new);
      }
      for(std::size_t me_adrs=0; me_adrs < E.size(); me_adrs++) {
	E[me_adrs] = tci::copy(ctx,E_new[me_adrs]);
      }
      RealT res_bp_err;
      BeliefPropagationCondition(ctx,edges_for_bp,
				 T,site_idx,site_to_mpi_rank,
				 E,edge_idx,comm,res_bp_err);
      if( res_bp_err < bp_tolerance ) {
	break;
      }
    }
    std::vector<BondDimT> res_bond_dim;
    std::vector<RealT> res_trunc_err;
    Truncation(ctx,edges_for_bp,T,site_idx,site_to_mpi_rank,
	       E,edge_idx,comm,max_bond_dim,sv_min,tg_err,
	       res_bond_dim,res_trunc_err);
  }
			  
  
  
}

#endif
