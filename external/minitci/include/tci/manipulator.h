/// This file is a part of r-ccs-cms/tnbp
/**
@file tci/manipulator.h
@brief header to declare manipulators
 */
#ifndef TCI_MANIPULATOR_H
#define TCI_MANIPULATOR_H
namespace tci {

  template <typename TenT>
  void set_elem(
       context_handle_t<TenT> &ctx,
       TenT &a,
       const elem_coors_t<TenT> &coors,
       const elem_t<TenT> el);

  template <typename TenT>
  void reshape(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const shape_t<TenT> &new_shape);

  template <typename TenT>
  void reshape(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const shape_t<TenT> &new_shape,
       TenT &out);

  template <typename TenT>
  void transpose(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const std::vector<bond_idx_t<TenT>> &new_order);

  template <typename TenT>
  void transpose(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const std::vector<bond_idx_t<TenT>> &new_order,
       TenT &out);

  template <typename TenT>
  void cplx_conj(
       context_handle_t<TenT> &ctx,
       TenT &inout);

  template <typename TenT>
  void cplx_conj(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       TenT &out);

  template <typename TenT>
  void to_cplx(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       cplx_ten_t<TenT> &out);

  template <typename TenT>
  cplx_ten_t<TenT> to_cplx(
       context_handle_t<TenT> &ctx,
       const TenT &in);

  template <typename TenT>
  void real(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       real_ten_t<TenT> &out);

  template <typename TenT>
  real_ten_t<TenT> real(
       context_handle_t<TenT> &ctx,
       const TenT &in);

  template <typename TenT>
  void imag(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       real_ten_t<TenT> &out);

  template <typename TenT>
  real_ten_t<TenT> imag(
       context_handle_t<TenT> &ctx,
       const TenT &in);

  template <typename TenT>
  void expand(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const std::map<bond_idx_t<TenT>, bond_dim_t<TenT>> &
       bond_idx_increment_map);

  template <typename TenT>
  void expand(
       context_handle_t<TenT> &ctx,
       TenT &in,
       const std::map<
       bond_idx_t<TenT>, bond_dim_t<TenT>> &
       bond_idx_increment_map,
       TenT &out);

  template <typename TenT>
  void shrink(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const Map<bond_idx_t<TenT>,Pair<elem_coor_t<TenT>,elem_coor_t<TenT>>> &
       bd_idx_el_coor_pair_map);

  template <typename TenT>
  void shrink(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const Map<bond_idx_t<TenT>,Pair<elem_coor_t<TenT>,elem_coor_t<TenT>>> &
       bd_idx_el_coor_pair_map,
       TenT &out);

  template <typename TenT>
  void extract_sub(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const std::vector<std::pair<
       elem_coor_t<TenT>,elem_coor_t<TenT>>> &
       coor_pairs);

  template <typename TenT>
  void extract_sub(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const std::vector<std::pair<
       elem_coor_t<TenT>,
       elem_coor_t<TenT>>> &
       coor_pairs,
       TenT &out);

  template <typename TenT>
  void replace_sub(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const TenT &sub,
       const elem_coors_t<TenT> &begin_pt);

  template <typename TenT>
  void replace_sub(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const TenT &sub,
       const elem_coors_t<TenT> &begin_pt,
       TenT &out);

  template <typename TenT>
  void concatenate(
       context_handle_t<TenT> &ctx,
       const std::vector<TenT> &ins,
       const bond_idx_t<TenT> concat_bdidx,
       TenT &out);

  template <typename TenT>
  void stack(
       context_handle_t<TenT> &ctx,
       const std::vector<TenT> &ins,
       const bond_idx_t<TenT> stack_bdidx,
       TenT &out);

  template <typename TenT, typename Func>
  void for_each(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       Func &&f);

  template <typename TenT, typename Func>
  TenT for_each(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       Func &&f);

  template <typename TenT, typename Func>
  void for_each_with_coors(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       Func &&f);

  template <typename TenT, typename Func>
  TenT for_each_with_coors(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       Func &&f);
  
}
#endif
