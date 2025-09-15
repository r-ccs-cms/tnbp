/// This file is a part of r-ccs-cms/tnbp
/**
@file tci/traits.h
@brief header to define traits
 */
#ifndef TCI_TRAITS_H
#define TCI_TRAITS_H
namespace tci {
  template <typename TenT>
  struct tensor_traits;

  template <typename TenT>
  using ten_t = typename tensor_traits<TenT>::ten_t;

  template <typename TenT>
  using rank_t = typename tensor_traits<TenT>::rank_t;

  template <typename TenT>
  using shape_t = typename tensor_traits<TenT>::shape_t;

  template <typename TenT>
  using bond_dim_t = typename tensor_traits<TenT>::bond_dim_t;

  template <typename TenT>
  using bond_idx_t = typename tensor_traits<TenT>::bond_idx_t;

  template <typename TenT>
  using bond_label_t = typename tensor_traits<TenT>::bond_label_t;

  template <typename TenT>
  using ten_size_t = typename tensor_traits<TenT>::ten_size_t;

  template <typename TenT>
  using elem_t = typename tensor_traits<TenT>::elem_t;

  template <typename TenT>
  using elem_coor_t = typename tensor_traits<TenT>::elem_coor_t;

  template <typename TenT>
  using elem_coors_t = typename tensor_traits<TenT>::elem_coors_t;

  template <typename TenT>
  using real_t = typename tensor_traits<TenT>::real_t;

  template <typename TenT>
  using real_ten_t = typename tensor_traits<TenT>::real_ten_t;

  template <typename TenT>
  using cplx_t = typename tensor_traits<TenT>::cplx_t;

  template <typename TenT>
  using cplx_ten_t = typename tensor_traits<TenT>::cplx_ten_t;

  template <typename TenT>
  using context_handle_t = typename tensor_traits<TenT>::context_handle_t;

  /**
     Alias declaration
   */
  template <typename ElemT>
  using List = std::vector<ElemT>;

  template <typename First, typename Second>
  using Pair = std::pair<First,Second>;

  template <typename Key, typename Value>
  using Map = std::map<Key,Value>;
  
    
}
#endif


