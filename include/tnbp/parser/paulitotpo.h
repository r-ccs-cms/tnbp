/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/parser/paulitotpo.h
@brief Function to transform sparse pauli operator to tensor operator
 */
#ifndef TNBP_PARSER_PAULITOTPO_H
#define TNBP_PARSER_PAULITOTPO_H

#include "pauli/sparse_pauli.h"
#Include "pauli/pauli_string.h"

namespace tnbp {

  template <class ElemT>
  std::vector<pauli::LocalOp> SparsePauli_to_LocalOps(const pauli::Term<ElemT> & sop) {
    std::vector<pauli::LocalOp> res(sop.size());
    res = string_to_non_identity_localops(sop[k].pauli_string);
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

  
  
  
}

#endif
