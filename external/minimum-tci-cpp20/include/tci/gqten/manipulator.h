/// This file is a part of r-ccs-cms/tnbp
/**
@file tci/gqten/manipulator.h
@brief header to define manipulator routines for gqten::tensor<ElemT>
*/
#ifndef TCI_GQTEN_MANIPULATOR_H
#define TCI_GQTEN_MANIPULATOR_H

#include <vector>
#include <map>
#include <complex>
#include <cstdlib>     // std::malloc, std::calloc
#include <utility>     // std::pair
#include <functional>  // std::invoke
#include <algorithm>   // std::min

namespace tci {

  /**
  template <typename TenT>
  void set_elem(
       context_handle_t<TenT> &ctx,
       TenT &a,
       const elem_coors_t<TenT> &coors,
       const elem_t<TenT> el);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void set_elem(
       context_handle_t<TenT> &ctx,
       TenT &a,
       const elem_coors_t<TenT> &coors,
       const elem_t<TenT> el) {
    a.SetElem(coors,el);
  }

  /**
  template <typename TenT>
  void reshape(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const shape_t<TenT> &new_shape);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void reshape(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const shape_t<TenT> &new_shape) {
    inout.Reshape(new_shape);
  }

  /**
  template <typename TenT>
  void reshape(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const shape_t<TenT> &new_shape,
       TenT &out);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void reshape(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const shape_t<TenT> &new_shape,
       TenT &out) {
    out = in;
    out.Reshape(new_shape);
  }

  /**
  template <typename TenT>
  void transpose(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const List<bond_idx_t<TenT>> &new_order);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void transpose(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const List<bond_idx_t<TenT>> &new_order) {
    inout.Transpose(new_order);
  }

  /**
  template <typename TenT>
  void transpose(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const List<bond_idx_t<TenT>> &new_order,
       TenT &out);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void transpose(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const List<bond_idx_t<TenT>> &new_order,
       TenT &out) {
    out = in;
    out.Transpose(new_order);
  }

  /**
  template <typename TenT>
  void cplx_conj(
       context_handle_t<TenT> &ctx,
       TenT &inout);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void cplx_conj(
       context_handle_t<TenT> &ctx,
       TenT &inout) {
    inout.Conjugate();
  }

  /**
  template <typename TenT>
  void cplx_conj(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       TenT &out);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void cplx_conj(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       TenT &out) {
    out = in;
    out.Conjugate();
  }

  /**
  template <typename TenT>
  void to_cplx(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       cplx_ten_t<TenT> &out);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void to_cplx(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       cplx_ten_t<TenT> &out) {
    using CplxTenT = typename tensor_traits<TenT>::cplx_ten_t;
    using CplxT = typename tensor_traits<TenT>::cplx_t;
    using ShapeT = typename tensor_traits<TenT>::shape_t;
    size_t size = in.Size();
    CplxT * praw = static_cast<CplxT*>(std::malloc(size*sizeof(CplxT)));
    const elem_t<TenT> * pvals = in.GetRaw();
    for(size_t k=0; k < size; k++) {
      praw[k] = static_cast<CplxT>(pvals[k]);
    }
    ShapeT shape = in.Shape();
    out = CplxTenT(shape,praw);
  }

  /**
  template <typename TenT>
  cplx_ten_t<TenT> to_cplx(
       context_handle_t<TenT> &ctx,
       const TenT &in);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  cplx_ten_t<TenT> to_cplx(
       context_handle_t<TenT> &ctx,
       const TenT &in) {
    using CplxTenT = typename tensor_traits<TenT>::cplx_ten_t;
    using CplxT = typename tensor_traits<TenT>::cplx_t;
    using ShapeT = typename tensor_traits<TenT>::shape_t;
    size_t size = in.Size();
    CplxT * praw = static_cast<CplxT*>(std::malloc(size*sizeof(CplxT)));
    const elem_t<TenT> * pvals = in.GetRaw();
    for(size_t k=0; k < size; k++) {
      praw[k] = static_cast<CplxT>(pvals[k]);
    }
    ShapeT shape = in.Shape();
    return CplxTenT(shape,praw);
  }

  /**
  template <typename TenT>
  void real(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       real_ten_t<TenT> &out);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void real(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       real_ten_t<TenT> &out) {
    using RealTenT = typename tensor_traits<TenT>::real_ten_t;
    using RealT = typename tensor_traits<TenT>::real_t;
    using ShapeT = typename tensor_traits<TenT>::shape_t;
    size_t size = in.Size();
    RealT * praw = static_cast<RealT*>(std::malloc(size*sizeof(RealT)));
    const elem_t<TenT> * pvals = in.GetRaw();
    for(size_t k=0; k < size; k++) {
      praw[k] = gqten::GetReal(pvals[k]);
    }
    ShapeT shape = in.Shape();
    out = RealTenT(shape,praw);
  }

  /**
  template <typename TenT>
  real_ten_t<TenT> real(
       context_handle_t<TenT> &ctx,
       const TenT &in);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  real_ten_t<TenT> real(
       context_handle_t<TenT> &ctx,
       const TenT &in) {
    using RealTenT = typename tensor_traits<TenT>::real_ten_t;
    using RealT = typename tensor_traits<TenT>::real_t;
    using ShapeT = typename tensor_traits<TenT>::shape_t;
    size_t size = in.Size();
    RealT * praw = static_cast<RealT*>(std::malloc(size*sizeof(RealT)));
    const elem_t<TenT> * pvals = in.GetRaw();
    for(size_t k=0; k < size; k++) {
      praw[k] = gqten::GetReal(pvals[k]);
    }
    ShapeT shape = in.Shape();
    return RealTenT(shape,praw);
  }


  inline float GetImag(const float a) { return 0.0; }
  inline double GetImag(const double a) { return 0.0; }
  inline float GetImag(const std::complex<float> a) { return a.imag(); }
  inline double GetImag(const std::complex<double> a) { return a.imag(); }
  /**
  template <typename TenT>
  void imag(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       real_ten_t<TenT> &out);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void imag(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       real_ten_t<TenT> &out) {
    using RealTenT = typename tensor_traits<TenT>::real_ten_t;
    using RealT = typename tensor_traits<TenT>::real_t;
    using ShapeT = typename tensor_traits<TenT>::shape_t;
    size_t size = in.Size();
    RealT * praw = static_cast<RealT*>(std::malloc(size*sizeof(RealT)));
    const elem_t<TenT> * pvals = in.GetRaw();
    for(size_t k=0; k < size; k++) {
      praw[k] = GetImag(pvals[k]);
    }
    ShapeT shape = in.Shape();
    out = RealTenT(shape,praw);
    
  }

  /**
  template <typename TenT>
  real_ten_t<TenT> imag(
       context_handle_t<TenT> &ctx,
       const TenT &in);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  real_ten_t<TenT> imag(
	context_handle_t<TenT> &ctx,
	const TenT &in) {
    using RealTenT = typename tensor_traits<TenT>::real_ten_t;
    using RealT = typename tensor_traits<TenT>::real_t;
    using ShapeT = typename tensor_traits<TenT>::shape_t;
    size_t size = in.Size();
    RealT * praw = static_cast<RealT*>(std::malloc(size*sizeof(RealT)));
    const elem_t<TenT> * pvals = in.GetRaw();
    for(size_t k=0; k < size; k++) {
      praw[k] = GetImag(pvals[k]);
    }
    ShapeT shape = in.Shape();
    return RealTenT(shape,praw);
  }

  /**
  template <typename TenT>
  void expand(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const Map<bond_idx_t<TenT>, bond_dim_t<TenT>> &
       bond_idx_increment_map);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void expand(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const Map<bond_idx_t<TenT>,
       bond_dim_t<TenT>> &
       bond_idx_increment_map) {
    using BondDimT = typename tensor_traits<TenT>::bond_dim_t;
    using BondIdxT = typename tensor_traits<TenT>::bond_idx_t;
    using ShapeT = typename tensor_traits<TenT>::shape_t;
    size_t size = inout.Size();
    ShapeT shape = inout.Shape();
    ShapeT new_shape = inout.Shape();
    
    for(const auto & [k,d] : bond_idx_increment_map) {
      new_shape[k] = shape[k] + d;
    }
    size_t new_size = 1;
    for(const auto & d : new_shape) {
      new_size *= d;
    }
    elem_t<TenT> * new_praw = static_cast<elem_t<TenT>*>(std::calloc(new_size,sizeof(elem_t<TenT>)));
    TenT temp(new_shape,new_praw);
    std::vector<BondIdxT> coordinate(inout.Rank());
    for(size_t i=0; i < size; i++) {
      size_t itemp = i;
      for(size_t l=0; l < inout.Rank(); l++) {
	coordinate[l] = itemp % shape[l];
	itemp /= shape[l];
      }
      temp.SetElem(coordinate,inout.GetElem(coordinate));
    }
    inout = std::move(temp);
  }

  /**
  template <typename TenT>
  void expand(
       context_handle_t<TenT> &ctx,
       TenT &in,
       const Map<
       bond_idx_t<TenT>, bond_dim_t<TenT>> &
       bond_idx_increment_map,
       TenT &out);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void expand(
       context_handle_t<TenT> &ctx,
       TenT &in,
       const Map<
       bond_idx_t<TenT>,
       bond_dim_t<TenT>> &
       bond_idx_increment_map,
       TenT &out) {
    using BondDimT = typename tensor_traits<TenT>::bond_dim_t;
    using BondIdxT = typename tensor_traits<TenT>::bond_idx_t;
    using ShapeT = typename tensor_traits<TenT>::shape_t;
    size_t size = in.Size();
    ShapeT shape = in.Shape();
    ShapeT new_shape = in.Shape();
    for(const auto & [k,d] : bond_idx_increment_map) {
      new_shape[k] = d;
    }
    size_t new_size = 1;
    for(const auto & d : new_shape) {
      new_size *= d;
    }
    elem_t<TenT> * new_praw = static_cast<elem_t<TenT>*>(std::calloc(new_size,sizeof(elem_t<TenT>)));
    out = TenT(new_shape,new_praw);
    std::vector<BondIdxT> coordinate(in.Rank());
    for(size_t i=0; i < size; i++) {
      size_t itemp = i;
      for(size_t l=0; l < in.Rank(); l++) {
	coordinate[l] = itemp % shape[l];
	itemp /= shape[l];
      }
      out.SetElem(coordinate,in.GetElem(coordinate));
    }
  }

  /**
  template <typename TenT>
  void shrink(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const bond_idx_elem_coor_pair_map<TenT> &
       bd_idx_el_coor_pair_map);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void shrink(
       context_handle_t<TenT>& /*ctx*/,
       TenT& inout,
       const Map<
       bond_idx_t<TenT>,
       Pair<elem_coor_t<TenT>,
            elem_coor_t<TenT>>>&
       bd_idx_el_coor_pair_map) {
    using Traits   = tensor_traits<TenT>;
    using BondIdxT = typename Traits::bond_idx_t;   // size_t
    using ShapeT   = typename Traits::shape_t;
    
    const ShapeT shape = inout.Shape();
    const size_t rank  = shape.size();
    
    ShapeT new_shape = shape;
    for (const auto& [axis, range] : bd_idx_el_coor_pair_map) {
      const BondIdxT a = static_cast<BondIdxT>(axis);
      if (a >= rank) throw std::out_of_range("shrink: axis out of range");
      const auto [first, second] = range;
      if (first > second || second > shape[a])
	throw std::out_of_range("shrink: invalid [first,second) for axis");
      new_shape[a] = static_cast<BondIdxT>(second - first);
    }
    
    size_t new_size = 1;
    for (auto d : new_shape) new_size *= d;
    elem_t<TenT>* raw = static_cast<elem_t<TenT>*>(std::calloc(new_size, sizeof(elem_t<TenT>)));
    if (!raw) throw std::bad_alloc();
    TenT out(new_shape, raw);
    
    std::vector<BondIdxT> src(rank), dst(rank);
    for (size_t i = 0; i < new_size; ++i) {
      size_t t = i;
      for (size_t ax = 0; ax < rank; ++ax) {
	dst[ax] = static_cast<BondIdxT>(t % new_shape[ax]);
	t      /= new_shape[ax];
      }
      for (size_t ax = 0; ax < rank; ++ax) {
	if (auto it = bd_idx_el_coor_pair_map.find(ax); it != bd_idx_el_coor_pair_map.end()) {
	  src[ax] = static_cast<BondIdxT>(dst[ax] + it->second.first);
	} else {
	  src[ax] = dst[ax];
	}
      }
      out.SetElem(dst, inout.GetElem(src));
    }
    
    inout = std::move(out);
  }

  /**
  template <typename TenT>
  void shrink(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const bond_idx_elem_coor_pair_map<TenT> &
       bd_idx_el_coor_pair_map,
       TenT &out);
  */
  
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void shrink(
       context_handle_t<TenT>& /*ctx*/,
       const TenT& in,
       const Map<
       bond_idx_t<TenT>,
       Pair<elem_coor_t<TenT>,
            elem_coor_t<TenT>>>&
       bd_idx_el_coor_pair_map,
       TenT & out) {
    out = in;
    shrink(*(context_handle_t<TenT>*)nullptr, out, bd_idx_el_coor_pair_map);
  }

  /**
  template <typename TenT>
  void extract_sub(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const List<Pair<
       elem_coor_t<TenT>,elem_coor_t<TenT>>> &
       coor_pairs);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void extract_sub(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const List<std::pair<
       elem_coor_t<TenT>,
       elem_coor_t<TenT>>> &
       coor_pairs) {

    using BondIdxT = bond_idx_t<TenT>;

    Map<BondIdxT, Pair<elem_coor_t<TenT>, elem_coor_t<TenT>>> m;
    m.clear();
    for (BondIdxT ax = 0; ax < coor_pairs.size(); ++ax) {
      const auto& [lo, hi] = coor_pairs[ax];
      if (lo != 0 || hi != inout.Shape()[ax]) {
	m.emplace(ax, std::make_pair(lo, hi));
      }
    }
    shrink(ctx, inout, m);
  }
  
  /**
  template <typename TenT>
  void extract_sub(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const List<Pair<
       elem_coor_t<TenT>,
       elem_coor_t<TenT>>> &
       coor_pairs,
       TenT &out);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void extract_sub(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const std::vector<Pair<
       elem_coor_t<TenT>,
       elem_coor_t<TenT>>> &
       coor_pairs,
       TenT &out) {
    out = in;
    extract_sub(ctx, out, coor_pairs);
  }

  /**
  template <typename TenT>
  void replace_sub(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const TenT &sub,
       const elem_coors_t<TenT> &begin_pt);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void replace_sub(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const TenT &sub,
       const elem_coors_t<TenT> &begin_pt) {
    using BondIdxT = typename tensor_traits<TenT>::bond_idx_t;
    using ShapeT = typename tensor_traits<TenT>::shape_t;
    size_t size = inout.Size();
    size_t sub_size = sub.Size();
    ShapeT shape = inout.Shape();
    ShapeT sub_shape = sub.Shape();
    std::vector<BondIdxT> coordinate(inout.Rank());
    std::vector<BondIdxT> sub_coordinate(sub.Rank());
    for(size_t i=0; i < sub_size; i++) {
      size_t itemp = i;
      for(size_t l=0; l < sub.Rank(); l++) {
	sub_coordinate[l] = itemp % sub_shape[l];
	coordinate[l] = sub_coordinate[l] + begin_pt[l];
	itemp /= sub_shape[l];
      }
      inout.SetElem(coordinate,sub.GetElem(sub_coordinate));
    }
  }

  /**
  template <typename TenT>
  void replace_sub(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const TenT &sub,
       const elem_coors_t<TenT> &begin_pt,
       TenT &out);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void replace_sub(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const TenT &sub,
       const elem_coors_t<TenT> &begin_pt,
       TenT &out) {
    using BondIdxT = typename tensor_traits<TenT>::bond_idx_t;
    using ShapeT = typename tensor_traits<TenT>::shape_t;
    out = in;
    size_t size = in.Size();
    size_t sub_size = sub.Size();
    ShapeT shape = in.Shape();
    ShapeT sub_shape = sub.Shape();
    std::vector<BondIdxT> coordinate(in.Rank());
    std::vector<BondIdxT> sub_coordinate(sub.Rank());
    for(size_t i=0; i < sub_size; i++) {
      size_t itemp = i;
      for(size_t l=0; l < sub.Rank(); l++) {
	sub_coordinate[l] = itemp % sub_shape[l];
	coordinate[l] = sub_coordinate[l] + begin_pt[l];
	itemp /= sub_shape[l];
      }
      out.SetElem(coordinate,sub.GetElem(sub_coordinate));
    }
  }

  /**
  template <typename TenT>
  void concatenate(
       context_handle_t<TenT> &ctx,
       const List<TenT> &ins,
       const bond_idx_t<TenT> concat_bdidx,
       TenT &out);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void concatenate(
       context_handle_t<TenT> &ctx,
       const std::vector<TenT> &ins,
       const bond_idx_t<TenT> concat_bdidx,
       TenT &out) {

    using BondIdxT = typename tensor_traits<TenT>::bond_idx_t;
    using ShapeT = typename tensor_traits<TenT>::shape_t;
    if( ins.size() == 0 ) {
      return;
    }
    std::vector<BondIdxT> init_pt(ins.size());
    ShapeT new_shape = ins[0].Shape();
    size_t new_size = 1;
    for(size_t k=0; k < new_shape.size(); k++) {
      if( k != concat_bdidx ) {
	new_size *= new_shape[k];
      }
    }
    init_pt[0] = 0;
    for(size_t m=1; m < init_pt.size(); m++) {
      ShapeT shape = ins[m-1].Shape();
      init_pt[m] = shape[concat_bdidx] + init_pt[m-1];
    }
    size_t new_cdim = 0;
    std::vector<ShapeT> shapes(ins.size());
    for(size_t m=0; m < init_pt.size(); m++) {
      shapes[m] = ins[m].Shape();
      if( m == 0 ) {
	init_pt[0] = 0;
      } else {
	init_pt[m] = shapes[m-1][concat_bdidx] + init_pt[m-1];
      }
      new_cdim += shapes[m][concat_bdidx];
    }
    new_shape[concat_bdidx] = new_cdim;
    new_size *= new_cdim;

    elem_t<TenT> * new_praw = static_cast<elem_t<TenT>*>(std::malloc(new_size*sizeof(elem_t<TenT>)));
    out = TenT(new_shape,new_praw);
    std::vector<BondIdxT> new_coordinate(out.Rank());
    std::vector<BondIdxT> coordinate(out.Rank());
    for(size_t m=0; m < init_pt.size(); m++) {
      for(size_t i=0; i < ins[m].Size(); i++) {
	size_t itemp = i;
	for(size_t l=0; l < out.Rank(); l++) {
	  coordinate[l] = itemp % shapes[m][l];
	  if( l == concat_bdidx ) {
	    new_coordinate[l] = coordinate[l] + init_pt[m];
	  } else {
	    new_coordinate[l] = coordinate[l];
	  }
	  itemp /= shapes[m][l];
	}
	out.SetElem(new_coordinate,ins[m].GetElem(coordinate));
      }
    }
    
  }

  /**
  template <typename TenT>
  void stack(
       context_handle_t<TenT> &ctx,
       const List<TenT> &ins,
       const bond_idx_t<TenT> stack_bdidx,
       TenT &out);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void stack(
       context_handle_t<TenT> &ctx,
       const std::vector<TenT> &ins,
       const bond_idx_t<TenT> stack_bdidx,
       TenT &out) {
    
    using BondIdxT = typename tensor_traits<TenT>::bond_idx_t;
    using ShapeT = typename tensor_traits<TenT>::shape_t;
    if( ins.size() == 0 ) {
      return;
    }
    ShapeT shape = ins[0].Shape();
    ShapeT new_shape(ins[0].Rank()+1);
    for (size_t l = 0; l < new_shape.size(); ++l) {
      if (l < stack_bdidx) new_shape[l] = shape[l];
      else if (l == stack_bdidx) new_shape[l] = static_cast<BondIdxT>(ins.size());
      else new_shape[l] = shape[l - 1];
    }
    size_t size = ins[0].Size();
    size_t new_size = ins.size() * size;
    elem_t<TenT> * new_praw = static_cast<elem_t<TenT>*>(std::malloc(new_size*sizeof(elem_t<TenT>)));
    out = TenT(new_shape,new_praw);
    std::vector<BondIdxT> new_coordinate(out.Rank());
    std::vector<BondIdxT> coordinate(ins[0].Rank());
    for(size_t m=0; m < ins.size(); m++) {
      for(size_t i=0; i < ins[m].Size(); i++) {
	size_t itemp = i;
	for(size_t l=0; l < out.Rank(); l++) {
	  if( l < stack_bdidx ) {
	    coordinate[l] = itemp % shape[l];
	    new_coordinate[l] = coordinate[l];
	    itemp /= shape[l];
	  } else if ( l > stack_bdidx ) {
	    coordinate[l-1] = itemp % shape[l-1];
	    new_coordinate[l] = coordinate[l];
	    itemp /= shape[l];
	  } else {
	    new_coordinate[l] = m;
	  }
	}
	out.SetElem(new_coordinate,ins[m].GetElem(coordinate));
      }
    }
  }

  /**
  template <typename TenT, typename Func>
  void for_each(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       Func &&f);
  */
  template <typename TenT, typename Func>
  requires is_gqten_tensor_v<TenT>
  void for_each(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       Func &&f) {
    for(size_t i=0; i < inout.Size(); i++) {
      elem_t<TenT> elem = inout.GetElem(i);
      std::invoke(f, elem);
      inout.SetElem(i,elem);
    }
  }

  /**
  template <typename TenT, typename Func>
  TenT for_each(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       Func &&f);
  */
  template <typename TenT, typename Func>
  requires is_gqten_tensor_v<TenT>
  TenT for_each(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       Func &&f) {
    TenT res(in);
    for(size_t i=0; i < res.Size(); i++) {
      elem_t<TenT> elem = res.GetElem(i);
      std::invoke(f, elem);
      res.SetElem(i,elem);
    }
    return res;
  }

  /**
  template <typename TenT, typename Func>
  void for_each_with_coors(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       Func &&f);
  */
  template <typename TenT, typename Func>
  requires is_gqten_tensor_v<TenT>
  void for_each_with_coors(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       Func &&f) {
    using ShapeT = typename tensor_traits<TenT>::shape_t;
    using CoorsT = typename tensor_traits<TenT>::elem_coors_t;
    CoorsT coors(inout.Rank());
    ShapeT shape = inout.Shape();
    for(size_t i=0; i < inout.Size(); i++) {
      size_t itemp = i;
      for(size_t l=0; l < inout.Rank(); l++) {
	coors[l] = itemp % shape[l];
	itemp /= shape[l];
      }
      auto elem = inout.GetElem(coors);
      std::invoke(f,elem,coors);
      inout.SetElem(coors,elem);
    }
  }

  /**
  template <typename TenT, typename Func>
  TenT for_each_with_coors(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       Func &&f);
  */
  template <typename TenT, typename Func>
  requires is_gqten_tensor_v<TenT>
  TenT  for_each_with_coors(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       Func &&f) {
    using ShapeT = typename tensor_traits<TenT>::shape_t;
    using CoorsT = typename tensor_traits<TenT>::elem_coors_t;
    TenT res(in);
    CoorsT coors(res.Rank());
    ShapeT shape = res.Shape();
    for(size_t i=0; i < res.Size(); i++) {
      size_t itemp = i;
      for(size_t l=0; l < res.Rank(); l++) {
	coors[l] = itemp % shape[l];
	itemp /= shape[l];
      }
      auto elem = in.GetElem(coors);
      std::invoke(f,elem,coors);
      res.SetElem(coors,elem);
    }
    return res;
  }

}

#endif
