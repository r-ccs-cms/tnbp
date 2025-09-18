/// This file is a part of r-ccs-cms/tnbp
/**
@file tci/linalg.h
@brief header to declare linear algebra functions
 */
#ifndef TCI_LINALG_H
#define TCI_LINALG_H
namespace tci {

  template <typename TenT>
  void diag(
       context_handle_t<TenT> &ctx,
       TenT &inout);

  template <typename TenT>
  void diag(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       TenT &out);

  template <typename TenT>
  real_t<TenT> norm(
       context_handle_t<TenT> &ctx,
       const TenT &a);

  template <typename TenT>
  real_t<TenT> normalize(
       context_handle_t<TenT> &ctx,
       TenT &inout);

  template <typename TenT>
  real_t<TenT> normalize(
       context_handle_t<TenT> &ctx,
       const TenT in,
       TenT &out);

  template <typename TenT>
  void scale(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const elem_t<TenT> s);

  template <typename TenT>
  void scale(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const elem_t<TenT> s,
       TenT &out);

  template <typename TenT>
  void trace(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const List<Pair<
       bond_idx_t<TenT>,
       bond_idx_t<TenT>>>
       &bdidx_pairs);

  template <typename TenT>
  void trace(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const List<Pair<
       bond_idx_t<TenT>,
       bond_idx_t<TenT>>>
       &bdidx_pairs,
       TenT &out);

  template <typename TenT>
  void exp(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const rank_t<TenT> num_of_bonds_as_rows);

  template <typename TenT>
  void exp(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const rank_t<TenT> num_of_bonds_as_rows,
       TenT &out);

  template <typename TenT>
  void inverse(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const rank_t<TenT> num_of_bonds_as_rows);

  template <typename TenT>
  void inverse(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const rank_t<TenT> num_of_bonds_as_rows,
       TenT &out);

  template <typename TenT>
  void contract(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const std::vector<bond_label_t<TenT>> &bd_labs_a,
       const TenT &b,
       const std::vector<bond_label_t<TenT>> &bd_labs_b,
       TenT &c,
       const std::vector<bond_label_t<TenT>> &bd_labs_c);

  template <typename TenT>
  void contract(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const std::string_view bd_labs_str_a,
       const TenT &b,
       const std::string_view bd_labs_str_b,
       TenT &c,
       const std::string_view bd_labs_str_c);

  template <typename TenT>
  void linear_combine(
       context_handle_t<TenT> &ctx,
       const std::vector<TenT> &ins,
       TenT &out);

  template <typename TenT>
  void linear_combine(
       context_handle_t<TenT> &ctx,
       const std::vector<TenT> &ins,
       const std::vector<elem_t<TenT>> &coefs,
       TenT &out);

  template <typename TenT>
  void svd(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const rank_t<TenT> &num_of_bonds_as_rows,
       TenT &u,
       real_ten_t<TenT> &s_diag,
       TenT &v_dag);

  template <typename TenT>
  void trunc_svd(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const rank_t<TenT> &num_of_bonds_as_rows,
       TenT &u,
       real_ten_t<TenT> &s_diag,
       TenT &v_dag,
       real_t<TenT> &trunc_err,
       const real_t<TenT> s_min);

  template <typename TenT>
  void trunc_svd(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const rank_t<TenT> &num_of_bonds_as_rows,
       TenT &u,
       real_ten_t<TenT> &s_diag,
       TenT &v_dag,
       real_t<TenT> &trunc_err,
       const bond_dim_t<TenT> chi_max,
       const real_t<TenT> s_min);

  template <typename TenT>
  void trunc_svd(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const rank_t<TenT> &num_of_bonds_as_rows, TenT &u,
       real_ten_t<TenT> &s_diag,
       TenT &v_dag,
       real_t<TenT> &trunc_err,
       const bond_dim_t<TenT> chi_min,
       const bond_dim_t<TenT> chi_max,
       const real_t<TenT> target_trunc_err,
       const real_t<TenT> s_min);

  template <typename TenT>
  void qr(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const rank_t<TenT> &num_of_bonds_as_rows,
       TenT &q,
       TenT &r);

  template <typename TenT>
  void lq(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const rank_t<TenT> &num_of_bonds_as_rows,
       TenT &l,
       TenT &q);

  template <typename TenT>
  void eigvals(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const rank_t<TenT> &num_of_bonds_as_rows,
       cplx_ten_t<TenT> &w_diag);

  template <typename TenT>
  void eigvalsh(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const rank_t<TenT> &num_of_bonds_as_rows,
       real_ten_t<TenT> &w_diag);

  template <typename TenT>
  void eig(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const rank_t<TenT> &num_of_bonds_as_rows,
       cplx_ten_t<TenT> &w_diag,
       cplx_ten_t<TenT> &v);

  template <typename TenT>
  void eigh(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const rank_t<TenT> &num_of_bonds_as_rows,
       real_ten_t<TenT> &w_diag,
       TenT &v);
  
}

#endif
