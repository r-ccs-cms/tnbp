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
    auto edges_orig = edges;
#ifdef BMPS_DEBUG
    if( mpi_rank == 0 ) {
      std::cout << " edges in bmps_increment_line before contraction:";
      for(const auto & [u,v] : edges) {
	  std::cout << " (" << u << "," << v << ")";
      }
      std::cout << std::endl;
    }
#endif
    for(std::size_t adrs_line_in = min_adrs_line_in;
	adrs_line_in <= max_adrs_line_in;
	++adrs_line_in) {
      auto site_in = line_in[adrs_line_in];
      auto bond_orig = GetSurroundingBondIndex(site_in,edges_orig);
      std::vector<IntT> bond_otype(bond_orig.size(),2);
      for(std::size_t k=0; k < bond_orig.size(); k++) {
	auto site_a = edges_orig[bond_orig[k]].first;
	auto site_b = edges_orig[bond_orig[k]].second;
	auto site_next = (site_a == site_in) ? site_b : site_a;
	auto it_is_zero = std::find(line_out.begin(),
					line_out.end(),
					site_next);
	if( it_is_zero != line_out.end() ) {
	  bond_otype[k] = 0;
	}
      }
      auto it_otype_zero = std::find(bond_otype.begin(),
				     bond_otype.end(),
				     0);
      auto bond_in = GetSurroundingBondIndex(site_in,edges);
      IntT target_edge_adrs;
      if( it_otype_zero != bond_otype.end() ) {
	auto k_orig = std::distance(bond_otype.begin(),
				    it_otype_zero);
	auto site_orig_a = edges_orig[bond_orig[k_orig]].first;
	auto site_orig_b = edges_orig[bond_orig[k_orig]].second;
	auto site_orig_next = (site_orig_a == site_in) ? site_orig_b : site_orig_a;
	std::size_t target_bond_index;
	for(std::size_t k=0; k < bond_in.size(); k++) {
	  auto site_a = edges[bond_in[k]].first;
	  auto site_b = edges[bond_in[k]].second;
	  auto site_next = (site_a == site_in) ? site_b : site_a;
	  if( site_next == site_orig_next ) {
	    target_bond_index = k;
	    break;
	  }
	}
	target_edge_adrs = bond_in[target_bond_index];
      } else {
	std::vector<IntT> bond_type(bond_in.size(),2);
	for(std::size_t k=0; k < bond_in.size(); k++) {
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
	if( it_type_zero != bond_type.end() ) {
	  auto target_bond_index = std::distance(bond_type.begin(),
						 it_type_zero);
	  target_edge_adrs = bond_in[target_bond_index];
	} else if (it_type_one != bond_type.end()) {
	  auto target_bond_index = std::distance(bond_type.begin(),
						 it_type_one);
	  target_edge_adrs = bond_in[target_bond_index];
	}
      }
      auto site_a = edges[target_edge_adrs].first;
      auto site_b = edges[target_edge_adrs].second;
      auto site_next = (site_a == site_in) ? site_b : site_a;
#ifdef BMPS_DEBUG
      if( mpi_rank == 0 ) {
	std::cout << " start contraction between "
		  << site_in << " site and "
		  << site_next << " site:";
      }
      auto mpi_rank_in = site_to_mpi_rank.at(site_in);
      auto mpi_rank_next = site_to_mpi_rank.at(site_next);
      if( mpi_rank == mpi_rank_in ) {
	auto it_site_idx_adrs = std::find(site_idx.begin(),
					  site_idx.end(),
					  site_in);
	auto site_idx_adrs = std::distance(site_idx.begin(),
					   it_site_idx_adrs);
	auto shape_in = tcapi::shape(ctx,T[site_idx_adrs]);
	std::cout << " tensor at site " << site_in << ":";
	for(auto const & dim : shape_in) {
	  std::cout << " " << dim;
	}
	std::cout << std::endl;
      }
      if( mpi_rank == mpi_rank_next ) {
	auto it_site_idx_adrs = std::find(site_idx.begin(),
					  site_idx.end(),
					  site_next);
	auto site_idx_adrs = std::distance(site_idx.begin(),
					   it_site_idx_adrs);
	auto shape_next = tcapi::shape(ctx,T[site_idx_adrs]);
	std::cout << " tensor at site " << site_next << ":";
	for(auto const & dim : shape_next) {
	  std::cout << " " << dim;
	}
	std::cout << std::endl;
      }
#endif
      graph_tensor_contraction(ctx,edges,T,site_idx,site_to_mpi_rank,
			       site_next,site_in,site_next,comm);
#ifdef BMPS_DEBUG
      if( mpi_rank == 0 ) {
	std::cout << " edges in bmps_increment_line after contraction between "
		  << site_in << " site and " << site_next << " site:";
	for(const auto & [u,v] : edges) {
	  std::cout << " (" << u << "," << v << ")";
	}
	std::cout << std::endl;
      }
#endif
    }
  }
  
  template <typename TenT, typename IntT>
  void bmps_increment_line_for_edges(const std::vector<IntT> & line_in,
				     const std::vector<IntT> & line_out,
				     std::vector<std::pair<IntT,IntT>> & edges) {
    auto edges_orig = edges;
    for(size_t adrs_line_in=0;
	adrs_line_in < line_in.size();
	++adrs_line_in) {
      auto site_in = line_in[adrs_line_in];
      auto bond_orig = GetSurroundingBondIndex(site_in,edges_orig);
      std::vector<IntT> bond_otype(bond_orig.size(),2);
      for(std::size_t k=0; k < bond_orig.size(); k++) {
	auto site_a = edges_orig[bond_orig[k]].first;
	auto site_b = edges_orig[bond_orig[k]].second;
	auto site_next = (site_a == site_in) ? site_b : site_a;
	auto it_is_zero = std::find(line_out.begin(),
					line_out.end(),
					site_next);
	if( it_is_zero != line_out.end() ) {
	  bond_otype[k] = 0;
	}
      }
      auto it_otype_zero = std::find(bond_otype.begin(),
				     bond_otype.end(),
				     0);
      auto bond_in = GetSurroundingBondIndex(site_in,edges);
      IntT target_edge_adrs;
      if( it_otype_zero != bond_otype.end() ) {
	auto k_orig = std::distance(bond_otype.begin(),
				    it_otype_zero);
	auto site_orig_a = edges_orig[bond_orig[k_orig]].first;
	auto site_orig_b = edges_orig[bond_orig[k_orig]].second;
	auto site_orig_next = (site_orig_a == site_in) ? site_orig_b : site_orig_a;
	std::size_t target_bond_index;
	for(std::size_t k=0; k < bond_in.size(); k++) {
	  auto site_a = edges[bond_in[k]].first;
	  auto site_b = edges[bond_in[k]].second;
	  auto site_next = (site_a == site_in) ? site_b : site_a;
	  if( site_next == site_orig_next ) {
	    target_bond_index = k;
	    break;
	  }
	}
	target_edge_adrs = bond_in[target_bond_index];
      } else {
	std::vector<IntT> bond_type(bond_in.size(),2);
	for(std::size_t k=0; k < bond_in.size(); k++) {
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
	if( it_type_zero != bond_type.end() ) {
	  auto target_bond_index = std::distance(bond_type.begin(),
						 it_type_zero);
	  target_edge_adrs = bond_in[target_bond_index];
	} else if (it_type_one != bond_type.end()) {
	  auto target_bond_index = std::distance(bond_type.begin(),
						 it_type_one);
	  target_edge_adrs = bond_in[target_bond_index];
	}
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
	E[me_adrs] = tcapi::copy(ctx,E_new[me_adrs]);
      }
      RealT res_bp_err;
      BeliefPropagationCondition(ctx,edges_for_bp,
				 T,site_idx,site_to_mpi_rank,
				 E,edge_idx,comm,res_bp_err);
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
