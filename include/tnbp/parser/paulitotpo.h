/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/parser/paulitotpo.h
@brief Function to transform sparse pauli operator to tensor operator
 */
#ifndef TNBP_PARSER_PAULITOTPO_H
#define TNBP_PARSER_PAULITOTPO_H

#include "pauli/sparse_pauli.h"
#include "pauli/pauli_string.h"

namespace tnbp {

  template <typename ElemT>
  std::vector<pauli::LocalOp> PauliTerm_to_LocalOps(const pauli::Term<ElemT> & sop) {
    std::vector<pauli::LocalOp> res;
    res = pauli::string_to_non_identity_localops(sop.pauli_string);
    return res;
  }
  
  template <typename ElemT>
  std::vector<ElemT> PauliOpMatrix(const pauli::LocalOp & p) {
    switch (p.op) {
    case 'I':
      return { ElemT(1,0), ElemT(0,0), ElemT(0,0), ElemT(1,0) };
    case 'X':
      return { ElemT(0,0), ElemT(1,0), ElemT(1,0), ElemT(0,0) };
    case 'Y':
      return { ElemT(0,0), ElemT(0,-1), ElemT(0,1), ElemT(0,0) };
    case 'Z':
      return { ElemT(1,0), ElemT(0,0), ElemT(0,0), ElemT(-1,0) };
    default:
      throw std::runtime_error("matrix for this op not implemented");
    }
  }

  inline size_t PauliOpSite(const pauli::LocalOp & p) { return p.qubit; }

  template <typename TenT>
  void SparsePauliToTensorOp(
       context_handle_t<TenT> & ctx,
       const std::vector<pauli::Term<elem_t<TenT>>> & sparse_pauli,
       std::vector<TenT> & res_tensor,
       std::vector<std::vector<int>> & res_qubits) {
    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using RealT = typename tci::tensor_traits<TenT>::real_t;
    using BondLabelT = typename tci::tensor_traits<TenT>::bond_label_t;
    using BondIdxT = typename tci::tensor_traits<TenT>::bond_idx_t;
    using RankT = typename tci::tensor_traits<TenT>::rank_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    using CoorsT = typename tci::tensor_traits<TenT>::elem_coors_t;

    res_tensor.resize(sparse_pauli.size());
    res_qubits.resize(sparse_pauli.size());
    auto it_tensor = res_tensor.begin();
    auto it_qubits = res_qubits.begin();
    for(const auto & term : sparse_pauli) {
      std::vector<pauli::LocalOp> localop =
	PauliTerm_to_LocalOps<ElemT>(term);
      std::vector<int> qubits;
      TenT opten;
      TenT tempo;
      TenT local;
      for(size_t k=0; k < localop.size(); k++) {
	qubits.push_back(PauliOpSite(localop[k]));
	ShapeT shapeL(2,2);
	auto mat = PauliOpMatrix<ElemT>(localop[k]);
	auto itmat = mat.begin();
	tci::assign_from_container(
	     ctx,shapeL,itmat,
	     [&shapeL](const CoorsT & coors) {
	       return coors[0]+shapeL[0]*coors[1];
	     },local);
	RankT rank = tci::rank(ctx,opten);
	List<BondLabelT> labelO(2*(k+1));
	List<BondLabelT> labelL(2);
	List<BondLabelT> labelT(2*(k+2));
	BondLabelT label = 0;
	for(size_t s=0; s < k+1; s++) {
	  labelO[s] = label++;
	}
	labelL[0] = label++;
	for(size_t s=k+1; s < 2*(k+1); s++) {
	  labelO[s] = label++;
	}
	labelL[1] = label++;
	std::iota(labelT.begin(),labelT.end(),0);
	tci::contract(ctx,opten,labelO,local,labelL,tempo,labelT);
	tci::copy(ctx,tempo,opten);
      }
      *it_qubits++ = qubits;
      tci::copy(ctx,opten,*it_tensor);
      it_tensor++;
    }
  }
	     
  template <typename TenT>
  void ClassifyOps(context_handle_t<TenT> & ctx,
		   const std::vector<std::pair<int,int>> & edges,
		   const std::vector<TenT> & tensor,
		   const std::vector<std::vector<int>> & qubits,
		   std::vector<TenT> & onesite_tensor,
		   std::vector<int>  & onesite,
		   std::vector<TenT> & twosite_tensor,
		   std::vector<std::pair<int,int>> & twosite) {

    std::vector<size_t> onesite_address;
    std::vector<size_t> twosite_address;
    for(size_t k=0; k < qubits.size(); k++) {
      if ( qubits[k].size() == 1 ) {
	onesite_address.push_back(k);
      } else if ( qubits[k].size() == 2 ) {
	std::vector<int> bond = GetSurroundingBondIndex(qubits[k][0],edges);
	bool find_edge = false;
	for(size_t m=0; m < bond.size(); m++) {
	  if( edges[bond[m]].first == qubits[k][0]
	      && edges[bond[m]].second == qubits[k][1] ) {
	    find_edge = true;
	    break;
	  }
	  if( edges[bond[m]].first == qubits[k][1]
	      && edges[bond[m]].second == qubits[k][0] ) {
	    find_edge = true;
	    break;
	  }
	}
	if( find_edge ) {
	  twosite_address.push_back(k);
	}
      }
    }
    size_t num_onesite = onesite_address.size();
    size_t num_twosite = twosite_address.size();
    onesite_tensor.resize(num_onesite);
    onesite.resize(num_onesite);
    twosite_tensor.resize(num_twosite);
    twosite.resize(num_twosite);
    auto it_one_site = onesite.begin();
    auto it_one_tensor = onesite_tensor.begin();
    for(size_t i=0; i < onesite_address.size(); i++) {
      tci::copy(ctx,tensor[onesite_address[i]],*it_one_tensor++);
      *it_one_site++ = qubits[onesite_address[i]][0];
    }
    auto it_two_site = twosite.begin();
    auto it_two_tensor = twosite_tensor.begin();
    for(size_t i=0; i < twosite_address.size(); i++) {
      *it_two_site++ = std::make_pair(qubits[twosite_address[i]][0],
				      qubits[twosite_address[i]][1]);
      tci::copy(ctx,tensor[twosite_address[i]],
		*it_two_tensor++);
    }
  }
		    
  
}

#endif
