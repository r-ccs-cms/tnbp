/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/framework/typedef.h
@brief type definitions in namespace tnbp
*/
#ifndef TNBP_FRAMEWORK_TYPEDEF_H
#define TNBP_FRAMEWORK_TYPEDEF_H

#include "tci/tci.h"

#include <vector>
#include <utility>
#include <map>

namespace tnbp {

  template <typename TenT>
  using rank_t =
    typename tci::tensor_traits<TenT>::rank_t;

  template <typename TenT>
  using shape_t =
    typename tci::tensor_traits<TenT>::shape_t;

  template <typename TenT>
  using bond_dim_t =
    typename tci::tensor_traits<TenT>::bond_dim_t;

  template <typename TenT>
  using bond_idx_t =
    typename tci::tensor_traits<TenT>::bond_idx_t;

  template <typename TenT>
  using bond_label_t =
    typename tci::tensor_traits<TenT>::bond_label_t;

  template <typename TenT>
  using ten_size_t =
    typename tci::tensor_traits<TenT>::ten_size_t;

  template <typename TenT>
  using elem_t =
    typename tci::tensor_traits<TenT>::elem_t;

  template <typename TenT>
  using elem_coor_t =
    typename tci::tensor_traits<TenT>::elem_coor_t;

  template <typename TenT>
  using elem_coors_t =
    typename tci::tensor_traits<TenT>::elem_coors_t;

  template <typename TenT>
  using real_t =
    typename tci::tensor_traits<TenT>::real_t;

  template <typename TenT>
  using real_ten_t =
    typename tci::tensor_traits<TenT>::real_ten_t;

  template <typename TenT>
  using cplx_t =
    typename tci::tensor_traits<TenT>::cplx_t;

  template <typename TenT>
  using cplx_ten_t =
    typename tci::tensor_traits<TenT>::cplx_ten_t;

  template <typename TenT>
  using context_handle_t =
    typename tci::tensor_traits<TenT>::context_handle_t;

  
  template <typename ElemT>
  using List = tci::List<ElemT>;

  template <typename First, typename Second>
  using Pair = tci::Pair<First,Second>;

  template <typename Key, typename Value>
  using Map = tci::Map<Key,Value>;
  
}

#endif
