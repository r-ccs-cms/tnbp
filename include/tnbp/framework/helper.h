/// This file is a part of r-ccs-cms/tnbp
/**
 */
#ifndef TNBP_FRAMEWORK_HELPER_H
#define TNBP_FRAMEWORK_HELPER_H

#include <vector>
#include <algorithm>

namespace tnbp {

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
  void tensor_contraction_label_helper(const tci::List<BondLabelT> & label_a,
				       const tci::List<BondLabelT> & label_b,
				       tci::List<BondLabelT> & label_c) {
    tci::List<BondLabelT> a = label_a;
    tci::List<BondLabelT> b = label_b;
    std::sort(a.begin(),a.end());
    std::sort(b.begin(),b.end());

    label_c.clear();

    std::set_symmetric_difference(
	 a.begin(),a.end(),
	 b.begin(),b.end(),
	 std::back_inserter(label_c));
  }


  /**
   * @brief A function for performing tensor contraction between neighboring nodes on a graph in a
   *        distributed environment, with automatic updates to the graph structure and associated tensors.
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
   * @param[in] site_c label of the site after contraction. site_idx, site_to_mpi_rank are also updated.
   *
   */
  template <typename TenT, typename IntT>
  void graph_tensor_contraction(context_handle_t<TenT> & ctx,
				std::vector<std::pair<IntT,IntT>> & edges,
				std::vector<TenT> & W,
				std::vector<IntT> & site_idx,
				std::map<IntT,IntT> & site_to_mpi_rank,
				IntT site_a, IntT site_b, IntT site_c,
				MPI_Comm comm) {

    using BondLabelT = typename tci::tensor_traits<TenT>::bond_label_t;
    using BondIdxT = typename tci::tensor_traits<TenT>::bond_idx_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;

    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    
    auto bond_a = GetSurroundingBondIndex(site_a,edges);
    auto bond_b = GetSurroundingBondIndex(site_b,edges);
    tci::List<BondLabelT> label_a(bond_a.size());
    tci::List<BondLabelT> label_b(bond_b.size());
    for(size_t m=0; m < label_a.size(); m++) {
      label_a[m] = static_cast<BondLabelT>(bond_a[m]);
    }
    for(size_t m=0; m < label_b.size(); m++) {
      label_b[m] = static_cast<BondLabelT>(bond_b[m]);
    }
    tci::List<BondLabelT> label_c;
    tensor_contraction_label_helper(label_a,label_b,label_c);

    int mpi_rank_a = static_cast<int>(site_to_mpi_rank[site_a]);
    int mpi_rank_b = static_cast<int>(site_to_mpi_rank[site_b]);
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
	B = tci::copy(ctx,W[adrs_b]);
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
	tci::contract(ctx,W[adrs_a],label_a,B,label_b,W[adrs_a],label_c);
      }
    }

    // preparation of edge compatible modification
    auto edges_temp = edges;
    auto edge_target = make_edge(site_a,site_b);
    std::vector<IntT> bond_c(label_c.size());
    for(size_t m=0; m < label_c.size(); m++) {
      bond_c[m] = static_cast<IntT>(label_c[m]);
    }

    // graph compatible modification
    for(size_t m=0; m < bond_b.size(); m++) {
      if( !same_edge(edges_temp[bond_b[m]],edge_target) ) {
	IntT site_x = edges_temp[bond_b[m]].first;
	IntT site_y = edges_temp[bond_b[m]].second;
	IntT site_p = (site_x == site_b) ? site_y : site_x;
	auto it_ap_edge = find_edge(edges_temp,site_a,site_p);
	auto it_bp_edge = find_edge(edges_temp,site_b,site_p);
	if( it_ap_edge != edges_temp.end() ) {
	  auto mpi_rank_p = site_to_mpi_rank[site_p];
	  auto ap_edge_label = std::distance(edges_temp.begin(),it_ap_edge);
	  auto bp_edge_label = std::distance(edges_temp.begin(),it_bp_edge);
	  int mpi_merge = 0;
	  if( mpi_rank_a == mpi_rank ) {
	    mpi_merge += 2;
	  }
	  if( mpi_rank_p == mpi_rank ) {
	    mpi_merge += 1;
	  }
	  if( mpi_merge == 2 || mpi_merge == 3 ) {
	    auto it_adrs_c = std::find(site_idx.begin(),site_idx.end(),site_a);
	    auto adrs_c = std::distance(site_idx.begin(),it_adrs_c);
	    auto order_c = tci::order(ctx,W[adrs_c]);
	    auto shape_c = tci::shape(ctx,W[adrs_c]);
	    auto it_ap_bond_idx = std::find(bond_c.begin(),bond_c.end(),ap_edge_label);
	    auto it_bp_bond_idx = std::find(bond_c.begin(),bond_c.end(),bp_edge_label);
	    auto ap_bond_idx = std::distance(bond_c.begin(),it_ap_bond_idx);
	    auto bp_bond_idx = std::distance(bond_c.begin(),it_bp_bond_idx);
	    tci::List<BondIdxT> new_label_c(order_c);
	    ShapeT new_shape_c(order_c-1);
	    std::iota(new_label_c.begin(),new_label_c.end(),0);
	    new_label_c.insert(new_label_c.begin()+ap_bond_idx+1,bp_bond_idx);
	    if( bp_bond_idx > ap_bond_idx ) {
	      new_label_c.erase(new_label_c.begin()+bp_bond_idx+1);
	    } else {
	      new_label_c.erase(new_label_c.begin()+bp_bond_idx);
	    }
	    for(size_t k=0; k < order_c; k++) {
	      if( new_label_c[k] < ap_bond_idx ) {
		new_shape_c[k] = shape_c[new_label_c[k]];
	      } else if ( new_label_c[k] > ap_bond_idx ) {
		new_shape_c[k-1] = shape_c[new_label_c[k]];
	      } else if ( new_label_c[k] == ap_bond_idx ) {
		new_shape_c[k] = shape_c[ap_bond_idx] * shape_c[bp_bond_idx];
		k++;
	      }
	    }
	    tci::transpose(ctx,W[adrs_c],new_label_c);
	    tci::reshape(ctx,W[adrs_c],new_shape_c);
	    // remove bp_edge_label from bond_c
	    bond_c.erase(it_bp_bond_idx);
	  }
	  if( mpi_merge == 1 || mpi_merge == 3 ) {
	    auto it_adrs_p = std::find(site_idx.begin(),site_idx.end(),site_p);
	    auto adrs_p = std::distance(site_idx.begin(),it_adrs_p);
	    auto order_p = tci::order(ctx,W[adrs_p]);
	    auto shape_p = tci::shape(ctx,W[adrs_p]);
	    auto bond_p = GetSurroundingBondIndex(site_p,edges_temp);
	    auto it_ap_bond_idx = std::find(bond_p.begin(),bond_p.end(),ap_edge_label);
	    auto it_bp_bond_idx = std::find(bond_p.begin(),bond_p.end(),bp_edge_label);
	    auto ap_bond_idx = std::distance(bond_p.begin(),it_ap_bond_idx);
	    auto bp_bond_idx = std::distance(bond_p.begin(),it_bp_bond_idx);
	    tci::List<BondIdxT> new_label_p(order_p);
	    ShapeT new_shape_p(order_p-1);
	    std::iota(new_label_p.begin(),new_label_p.end(),0);
	    new_label_p.insert(new_label_p.begin()+ap_bond_idx+1,bp_bond_idx);
	    if( bp_bond_idx > ap_bond_idx ) {
	      new_label_p.erase(new_label_p.begin()+bp_bond_idx+1);
	    } else {
	      new_label_p.erase(new_label_p.begin()+bp_bond_idx);
	    }
	    for(size_t k=0; k < order_p; k++) {
	      if( new_label_p[k] < ap_bond_idx ) {
		new_shape_p[k] = shape_p[new_label_p[k]];
	      } else if ( new_label_p[k] > ap_bond_idx ) {
		new_shape_p[k-1] = shape_p[new_label_p[k]];
	      } else if ( new_label_p[k] == ap_bond_idx ) {
		new_shape_p[k] = shape_p[ap_bond_idx]*shape_p[bp_bond_idx];
		k++;
	      }
	    }
	    tci::transpose(ctx,W[adrs_p],new_label_p);
	    tci::reshape(ctx,W[adrs_p],new_shape_p);
	  }
	  // erase corresponding bp_bond_label
	  auto it_remove_edge = find_edge(edges,site_b,site_p);
	  edges.erase(it_remove_edge);
	  auto it_update_edge = find_edge(edges,site_a,site_p);
	  *it_update_edge = make_edge(site_a,site_p);
	} else if ( it_bp_edge != edges_temp.end() ) {
	  auto it_update_edge = find_edge(edges,site_b,site_p);
	  *it_update_edge = make_edge(site_a,site_p);
	}
      }
    }

    // up to here, site_c is treated as site_a
    auto it_target_edge = find_edge(edges,site_a,site_b);
    edges.erase(it_target_edge);
    // at here, site_b becomes absent in edges
    auto bond_update = GetSurroundingBondIndex(site_a,edges);
    for(size_t m=0; m < bond_update.size(); m++) {
      auto site_x = edges[bond_update[m]].first;
      auto site_y = edges[bond_update[m]].second;
      auto site_p = (site_x == site_a) ? site_y : site_x;
      edges[bond_update[m]] = make_edge(site_c,site_p);
    }

    auto it_site_a = std::find(site_idx.begin(),site_idx.end(),site_a);
    if( it_site_a != site_idx.end() ) {
      *it_site_a = site_c;
    }
    auto it_site_b = std::find(site_idx.begin(),site_idx.end(),site_b);
    if( it_site_b != site_idx.end() ) {
      auto adrs_b = std::distance(site_idx.begin(),it_site_b);
      auto it_w_b = W.begin();
      std::advance(it_w_b,adrs_b);
      W.erase(it_w_b);
      site_idx.erase(it_site_b);

    }

    auto node = site_to_mpi_rank.extract(site_a);
    if (!node.empty()) {
      node.key() = site_c;
      site_to_mpi_rank.insert(std::move(node));
    }
  }

  /**
   * @brief Edge update function same with graph_tensor_contraction without tensor contraction
   * @tparam IntT type of site label
   * @param[in/out] edges Edges to define the graph of current system
   * @param[in/out] site_a a label of a site which will be contracted
   * @param[in/out] site_b a label of another site which will be contracted
   * @param[in/out] site_c a label for updated tensor site
   *
   */

  template <typename IntT>
  void graph_tensor_contraction_for_edges(std::vector<std::pair<IntT,IntT>> & edges,
					  IntT site_a, IntT site_b, IntT site_c) {
    
    auto bond_a = GetSurroundingBondIndex(site_a,edges);
    auto bond_b = GetSurroundingBondIndex(site_b,edges);
    tci::List<IntT> label_a(bond_a.size());
    tci::List<IntT> label_b(bond_b.size());
    tci::List<IntT> label_c;
    for(size_t m=0; m < label_a.size(); m++) {
      label_a[m] = static_cast<IntT>(bond_a[m]);
    }
    for(size_t m=0; m < label_b.size(); m++) {
      label_b[m] = static_cast<IntT>(bond_b[m]);
    }
    tensor_contraction_label_helper(label_a,label_b,label_c);
    
    auto edges_temp = edges;
    auto edge_target = make_edge(site_a,site_b);
    std::vector<IntT> bond_c(label_c.size());
    for(size_t m=0; m < label_c.size(); m++) {
      bond_c[m] = static_cast<IntT>(label_c[m]);
    }
    // graph compatible modification
    for(size_t m=0; m < bond_b.size(); m++) {
      if( !same_edge(edges_temp[bond_b[m]],edge_target) ) {
	IntT site_x = edges_temp[bond_b[m]].first;
	IntT site_y = edges_temp[bond_b[m]].second;
	IntT site_p = (site_x == site_b) ? site_y : site_x;
	auto it_ap_edge = find_edge(edges_temp,site_a,site_p);
	auto it_bp_edge = find_edge(edges_temp,site_b,site_p);
	if( it_ap_edge != edges_temp.end() ) {
	  auto ap_edge_label = std::distance(edges_temp.begin(),it_ap_edge);
	  auto bp_edge_label = std::distance(edges_temp.begin(),it_bp_edge);
	  // erase corresponding bp_bond_label
	  auto it_remove_edge = find_edge(edges,site_b,site_p);
	  edges.erase(it_remove_edge);
	  auto it_update_edge = find_edge(edges,site_a,site_p);
	  *it_update_edge = make_edge(site_a,site_p);
	} else if ( it_bp_edge != edges_temp.end() ) {
	  auto it_update_edge = find_edge(edges,site_b,site_p);
	  *it_update_edge = make_edge(site_a,site_p);
	}
      }
    }
    auto it_target_edge = find_edge(edges,site_a,site_b);
    edges.erase(it_target_edge);
    auto bond_update = GetSurroundingBondIndex(site_a,edges);
    for(size_t m=0; m < bond_update.size(); m++) {
      auto site_x = edges[bond_update[m]].first;
      auto site_y = edges[bond_update[m]].second;
      auto site_p = (site_x == site_a) ? site_y : site_x;
      edges[bond_update[m]] = make_edge(site_c,site_p);
    }
  }
				
  
}

#endif
