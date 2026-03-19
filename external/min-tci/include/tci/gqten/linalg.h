#ifndef TCI_GQTEN_LINALG_H
#define TCI_GQTEN_LINALG_H

// C++17-compatible rewrite of linear algebra utilities for gqten backend.
// - replaces `requires is_gqten_tensor_v<TenT>` with SFINAE on template parameter
// - preserves original interfaces and semantics
// - depends on tci/traits.h providing: elem_t<T>, real_t<T>, cplx_t<T>, real_ten_t<T>, order_t<T>, bond_dim_t<T>,
//   bond_idx_t<T>, bond_label_t<T>, context_handle_t<T>, and is_gqten_tensor_v<T>.

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <algorithm>
#include <cstdlib>
#include <thread>
#include <chrono>
#include <string_view>
#include <cctype>
#include <string>
#include <limits>
#include <type_traits>
#include <numeric>   // std::iota
#include <cassert>   // assert

namespace tci {

  //----------------------------------------------------------------------------
  // diag (in-place)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void diag(
       context_handle_t<TenT> &ctx,
       TenT &inout) {
    size_t size = inout.Size();
    elem_t<TenT> * praw = static_cast<elem_t<TenT>*>(std::calloc(size*size,sizeof(elem_t<TenT>)));
    if(!praw) throw std::bad_alloc();
    const elem_t<TenT> * pdiag = inout.GetRaw();
    for(size_t i=0; i < size; i++) {
      praw[i+size*i] = pdiag[i];
    }
    inout = TenT({static_cast<int32_t>(size),
                  static_cast<int32_t>(size)},praw);
  }

  //----------------------------------------------------------------------------
  // diag (out-param)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void diag(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       TenT &out) {
    size_t dim = in.Size();
    elem_t<TenT> * praw = static_cast<elem_t<TenT>*>(std::calloc(dim*dim,sizeof(elem_t<TenT>)));
    if(!praw) throw std::bad_alloc();
    const elem_t<TenT> * pdiag = in.GetRaw();
    for(size_t i=0; i < dim; i++) {
      praw[i+dim*i] = pdiag[i];
    }
    out = TenT({dim,dim},praw);
  }

  //----------------------------------------------------------------------------
  // norm
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline real_t<TenT> norm(
       context_handle_t<TenT> &ctx,
       const TenT &a) {
    return a.CalcNorm();
  }

  //----------------------------------------------------------------------------
  // normalize (in-place)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline real_t<TenT> normalize(
       context_handle_t<TenT> &ctx,
       TenT &inout) {
    using RealT = typename tensor_traits<TenT>::real_t;
    RealT res = inout.Normalize();
    return res;
  }

  //----------------------------------------------------------------------------
  // normalize (out-param)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline real_t<TenT> normalize(
       context_handle_t<TenT> &ctx,
       const TenT in,
       TenT &out) {
    using RealT = typename tensor_traits<TenT>::real_t;
    out = in;
    RealT res = out.Normalize();
    return res;
  }

  //----------------------------------------------------------------------------
  // scale (in-place)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void scale(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const elem_t<TenT> s) {
    elem_t<TenT> value = inout.GetScale();
    inout.SetScale(value*s);
  }

  //----------------------------------------------------------------------------
  // scale (out-param)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void scale(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const elem_t<TenT> s,
       TenT &out) {
    out = in;
    elem_t<TenT> value = out.GetScale();
    out.SetScale(value*s);
  }

  //----------------------------------------------------------------------------
  // trace (in-place)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void trace(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const List<Pair<bond_idx_t<TenT>, bond_idx_t<TenT>>> & bdidx_pairs) {
    TenT temp;
    gqten::Trace(&inout,bdidx_pairs,&temp);
    inout = std::move(temp);
  }

  //----------------------------------------------------------------------------
  // trace (out-param)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void trace(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const List<Pair<bond_idx_t<TenT>, bond_idx_t<TenT>>> & bdidx_pairs,
       TenT & out) {
    gqten::Trace(&in,bdidx_pairs,&out);
  }

  //----------------------------------------------------------------------------
  // exp (in-place)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void exp(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const order_t<TenT> & num_of_bonds_as_rows) {
    TenT temp;
    gqten::ExpHermExact(&inout,num_of_bonds_as_rows,&temp);
    inout = std::move(temp);
  }

  //----------------------------------------------------------------------------
  // exp (out-param)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void exp(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const order_t<TenT> & num_of_bonds_as_rows,
       TenT &out) {
    gqten::ExpHermExact(&in,num_of_bonds_as_rows,&out);
  }

  //----------------------------------------------------------------------------
  // inverse (in-place)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void inverse(
       context_handle_t<TenT> &ctx,
       TenT &inout,
       const order_t<TenT> & num_of_bonds_as_rows) {
    TenT temp;
    gqten::Inverse(&inout,num_of_bonds_as_rows,&temp);
    inout = std::move(temp);
  }

  //----------------------------------------------------------------------------
  // inverse (out-param)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void inverse(
       context_handle_t<TenT> &ctx,
       const TenT &in,
       const order_t<TenT> & num_of_bonds_as_rows,
       TenT &out) {
    gqten::Inverse(&in,num_of_bonds_as_rows,&out);
  }

  //----------------------------------------------------------------------------
  // Label helpers
  //----------------------------------------------------------------------------
  template <typename BondLabelT>
  inline void GqtenCntLabelHelper(
    const List<BondLabelT> & bond_labs_a,
    const List<BondLabelT> & bond_labs_b,
    const List<BondLabelT> & bond_labs_c,
    std::vector<int> & labs_a,
    std::vector<int> & labs_b) {
    const auto nA = bond_labs_a.size();
    const auto nB = bond_labs_b.size();

    labs_a.resize(nA);
    labs_b.resize(nB);

    {
      std::unordered_set<int> seen;
      for (int x : bond_labs_a) {
        if (!seen.insert(x).second)
          throw std::invalid_argument("bond_labs_a has duplicate label.");
      }
    }
    {
      std::unordered_set<int> seen;
      for (int x : bond_labs_b) {
        if (!seen.insert(x).second)
          throw std::invalid_argument("bond_labs_b has duplicate label.");
      }
    }

    // --- Find common label between A and B ---
    std::unordered_set<int> setA(bond_labs_a.begin(), bond_labs_a.end());
    std::unordered_set<int> setB(bond_labs_b.begin(), bond_labs_b.end());
    std::vector<int> contracted;
    contracted.reserve(std::min(nA, nB));
    for (int x : setA) if (setB.count(x)) contracted.push_back(x);

    // --- find remaining labels for C
    //   survivorsExpected = (A or B) \ (A and  B)
    std::unordered_set<int> survivorsExpected = setA;
    survivorsExpected.insert(setB.begin(), setB.end());
    for (int x : contracted) survivorsExpected.erase(x);

    // assumption of bond_labs_c
    //  1) no same labels
    //  2) set should be equivalent to survivorsExpected
    {
      std::unordered_set<int> seen;
      for (int x : bond_labs_c) {
        if (!seen.insert(x).second)
          throw std::invalid_argument("bond_labs_c has duplicate label.");
      }
      if (seen != survivorsExpected)
        throw std::invalid_argument("bond_labs_c does not match the set of survivor labels.");
    }

    // --- assign labels in the order of bond_labs_c ---
    std::unordered_map<int,int> posMap; posMap.reserve(bond_labs_c.size());
    {
      int nextPos = 0;
      for (int x : bond_labs_c) posMap.emplace(x, nextPos++);
    }

    // --- assign negative labels for contraction ---
    //     order is defined in assending order
    std::sort(contracted.begin(), contracted.end());
    std::unordered_map<int,int> negMap; negMap.reserve(contracted.size());
    {
      int k = 1;
      for (int x : contracted) negMap.emplace(x, -k++);
    }

    // --- construct labs_a, labs_b ---
    auto convert_one = [&](const std::vector<int>& in,
                           std::vector<int>& out) {
      for (size_t i = 0; i < in.size(); ++i) {
        int lab = in[i];
        auto pit = posMap.find(lab);
        if (pit != posMap.end()) {      // remaining -> positive label
          out[i] = pit->second;
        } else {
          auto nit = negMap.find(lab);
          if (nit != negMap.end()) {    // contract -> negative label
            out[i] = nit->second;
          } else {
            // if remaining is absent in C, error
            throw std::logic_error("Label not found in pos/neg maps (inconsistent inputs).");
          }
        }
      }
    };

    convert_one(bond_labs_a, labs_a);
    convert_one(bond_labs_b, labs_b);

    // --- check no same labels in positive labels ---
    auto check_positive_unique = [](const std::vector<int>& v){
      std::unordered_set<int> pos;
      for (int x : v) if (x >= 0) if (!pos.insert(x).second)
        throw std::invalid_argument("Positive labels duplicated in output labs.");
    };
    check_positive_unique(labs_a);
    check_positive_unique(labs_b);
  }

  //----------------------------------------------------------------------------
  // contract (vector label version)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void contract(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const List<bond_label_t<TenT>> &bd_labs_a,
       const TenT &b,
       const List<bond_label_t<TenT>> &bd_labs_b,
       TenT &c,
       const List<bond_label_t<TenT>> &bd_labs_c) {
    std::vector<int> labs_a;
    std::vector<int> labs_b;
    GqtenCntLabelHelper(bd_labs_a,bd_labs_b,bd_labs_c,labs_a,labs_b);
    gqten::Contract(&a,&b,labs_a,labs_b,&c);
  }

  //----------------------------------------------------------------------------
  // Helpers for std::string_view
  //----------------------------------------------------------------------------
  template <typename BondLabelT>
  static inline void ParseBondLabelsOne(
        std::string_view s,
        std::unordered_map<std::string,BondLabelT>& sym2id, // common symbol table
        BondLabelT& nextId,                                  // next assign label
        std::vector<BondLabelT>& out) {
    auto is_sep = [](char c){
      return c==',' || std::isspace(static_cast<unsigned char>(c));
    };

    // Detect whether separator exists
    bool has_sep = false;
    for (char c : s) if (is_sep(c)) { has_sep = true; break; }

    auto push_token = [&](std::string tok){
      // neglect empty token
      if (tok.empty()) return;

      // allow to use it as it is if it is integer
      char* end=nullptr;
      long v = std::strtol(tok.c_str(), &end, 10);
      if (end && *end=='\0') {
        if (v < std::numeric_limits<BondLabelT>::min() ||
            v > std::numeric_limits<BondLabelT>::max())
          throw std::out_of_range("label integer out of int range");
        out.push_back(static_cast<int>(v));
        return;
      }
      // if it is not integer, assign integer using symbol table
      auto it = sym2id.find(tok);
      if (it == sym2id.end()) {
        sym2id.emplace(tok, nextId);
        out.push_back(nextId++);
      } else {
        out.push_back(it->second);
      }
    };

    if (has_sep) {
      // if there are separator: split token by comma/space
      std::string cur;
      cur.reserve(16);
      for (size_t i=0;i<s.size();++i) {
        char c = s[i];
        if (is_sep(c)) { push_token(cur); cur.clear(); }
        else           { cur.push_back(c); }
      }
      push_token(cur);
    } else {
      // if there are no separator: one symbol as one label (like "ij")
      // if you want to use "ab12" as label, it is better to use string with separators
      for (char c : s) {
        std::string tok(1, c);
        push_token(std::move(tok));
      }
    }
  }

  // Parse all labels for A/B/C at once, and obtain std::vector<int> using same map
  template <typename BondLabelT>
  static inline void ParseBondLabelsTriple(
         std::string_view a, std::string_view b, std::string_view c,
         List<BondLabelT>& out_a,
         List<BondLabelT>& out_b,
         List<BondLabelT>& out_c) {
    std::unordered_map<std::string,BondLabelT> sym2id;
    sym2id.reserve(32);
    BondLabelT nextId = 0;

    ParseBondLabelsOne(a, sym2id, nextId, out_a);
    ParseBondLabelsOne(b, sym2id, nextId, out_b);
    ParseBondLabelsOne(c, sym2id, nextId, out_c);
  }

  // Utility wrapper: string_view -> bond_labs -> labs for gqten
  inline void MakeGqtenContractLabelsFromStrings(
         std::string_view bd_labs_str_a,
         std::string_view bd_labs_str_b,
         std::string_view bd_labs_str_c,
         std::vector<int>& labs_a,     // output label for gqten::Contract
         std::vector<int>& labs_b) {   // output label for gqten::Contract
    std::vector<int> bond_labs_a, bond_labs_b, bond_labs_c;
    ParseBondLabelsTriple(bd_labs_str_a, bd_labs_str_b, bd_labs_str_c,
                          bond_labs_a, bond_labs_b, bond_labs_c);

    // Transform string label for A/B/C to (labs_a, labs_b) for gqten
    GqtenCntLabelHelper(bond_labs_a, bond_labs_b, bond_labs_c, labs_a, labs_b);
  }

  //----------------------------------------------------------------------------
  // contract (string label version)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void contract(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const std::string_view bd_labs_str_a,
       const TenT &b,
       const std::string_view bd_labs_str_b,
       TenT &c,
       const std::string_view bd_labs_str_c) {
    std::vector<int> labs_a;
    std::vector<int> labs_b;
    MakeGqtenContractLabelsFromStrings(bd_labs_str_a,bd_labs_str_b,bd_labs_str_c,
                                      labs_a,labs_b);
    gqten::Contract(&a,&b,labs_a,labs_b,&c);
  }

  //----------------------------------------------------------------------------
  // linear_combine (uniform coefficients)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void linear_combine(
       context_handle_t<TenT> & ctx,
       const List<TenT> & ins,
       TenT & out) {
    std::vector<TenT> ins_copy = ins;
    std::vector<TenT*> pins;
    pins.reserve(ins_copy.size());
    for (size_t i = 0; i < ins.size(); ++i) pins.emplace_back(&ins_copy[i]);
    std::vector<elem_t<TenT>> coef(ins.size(),static_cast<elem_t<TenT>>(1.0));
    gqten::LinearCombine(coef,pins,&out);
  }

  //----------------------------------------------------------------------------
  // linear_combine (uniform coefficients)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline TenT void linear_combine(
       context_handle_t<TenT> & ctx,
       const List<TenT> & ins) {
    TenT out;
    linear_combine(ctx,ins,out);
    return out;
  }

  //----------------------------------------------------------------------------
  // linear_combine (explicit coefficients)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void linear_combine(
       context_handle_t<TenT> &ctx,
       const List<TenT> &ins,
       const List<elem_t<TenT>> &coefs,
       TenT &out) {
    std::vector<TenT> ins_copy = ins;
    std::vector<TenT*> pins;
    pins.reserve(ins_copy.size());
    for (size_t i = 0; i < ins_copy.size(); ++i) {
      pins.emplace_back(&ins_copy[i]);
    }
    gqten::LinearCombine(coefs,pins,&out);
  }

  //----------------------------------------------------------------------------
  // linear_combine (explicit coefficients, return-value-type)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline TenT linear_combine(
       context_handle_t<TenT> &ctx,
       const List<TenT> &ins,
       const List<elem_t<TenT>> &coefs) {
    TenT out;
    linear_combine(ctx,ins,coefs,out);
    return out;
  }
  

  //----------------------------------------------------------------------------
  // svd
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void svd(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const order_t<TenT> &num_of_bonds_as_rows,
       TenT &u,
       real_ten_t<TenT> &s_diag,
       TenT &v_dag) {
    using RealT = typename tensor_traits<TenT>::real_t;
    RealT * ps_raw = nullptr;
    size_t k;
    gqten::SVD(&a,num_of_bonds_as_rows,&u,&v_dag,ps_raw,&k);
    s_diag = gqten::tensor<RealT>({static_cast<int32_t>(k)},ps_raw);
  }

  //----------------------------------------------------------------------------
  // trunc_svd (s_min only)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void trunc_svd(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const order_t<TenT> &num_of_bonds_as_rows,
       TenT &u,
       real_ten_t<TenT> &s_diag,
       TenT &v_dag,
       real_t<TenT> &trunc_err,
       const real_t<TenT> s_min) {
    using RealT = typename tensor_traits<TenT>::real_t;
    using ShapeT = typename tensor_traits<TenT>::shape_t;
    RealT * ps_raw = nullptr;
    size_t chi;
    const auto shape_a = a.Shape();
    const size_t order_a = a.Rank();
    size_t ldims = static_cast<size_t>(num_of_bonds_as_rows);

    size_t m = 1;
    size_t n = 1;
    for(size_t i=0; i < ldims; i++) m *= shape_a[i];
    for(size_t i=ldims; i < order_a; i++) n *= shape_a[i];

    size_t chi_max = std::min(m,n);

    gqten::TruncSVD(&a,ldims,chi_max,
                    &u,&v_dag,ps_raw,
                    &chi,&trunc_err,s_min);
    s_diag = gqten::tensor<RealT>({static_cast<int32_t>(chi)},ps_raw);
  }

  //----------------------------------------------------------------------------
  // trunc_svd (limit chi_max + s_min)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void trunc_svd(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const order_t<TenT> &num_of_bonds_as_rows,
       TenT &u,
       real_ten_t<TenT> &s_diag,
       TenT &v_dag,
       real_t<TenT> &trunc_err,
       const bond_dim_t<TenT> chi_max,
       const real_t<TenT> s_min) {
    using RealT = typename tensor_traits<TenT>::real_t;
    RealT * ps_raw;
    size_t chi;
    gqten::TruncSVD(&a,static_cast<size_t>(num_of_bonds_as_rows),
                    static_cast<size_t>(chi_max),
                    &u,&v_dag,ps_raw,&chi,&trunc_err,s_min);
    s_diag = gqten::tensor<RealT>({static_cast<int32_t>(chi)},ps_raw);
  }

  //----------------------------------------------------------------------------
  // trunc_svd (target error + chi range + s_min)
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void trunc_svd(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const order_t<TenT> &num_of_bonds_as_rows,
       TenT &u,
       real_ten_t<TenT> &s_diag,
       TenT &v_dag,
       real_t<TenT> &trunc_err,
       const bond_dim_t<TenT> chi_min,
       const bond_dim_t<TenT> chi_max,
       const real_t<TenT> target_trunc_err,
       const real_t<TenT> s_min) {
    using RealT = typename tensor_traits<TenT>::real_t;
    RealT * ps_raw;
    size_t chi;
    gqten::TruncSVD(&a,static_cast<size_t>(num_of_bonds_as_rows),
                    static_cast<RealT>(target_trunc_err),
                    static_cast<size_t>(chi_max),
                    static_cast<size_t>(chi_min),
                    &u,&v_dag,ps_raw,&chi,&trunc_err,s_min);
    s_diag = gqten::tensor<RealT>({static_cast<int32_t>(chi)},ps_raw);
  }

  //----------------------------------------------------------------------------
  // qr
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void qr(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const order_t<TenT> &num_of_bonds_as_rows,
       TenT &q,
       TenT &r) {
    gqten::QR(&a,num_of_bonds_as_rows,&q,&r);
  }

  //----------------------------------------------------------------------------
  // lq
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void lq(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const order_t<TenT> &num_of_bonds_as_rows,
       TenT &l,
       TenT &q) {
    using OrderT = typename tensor_traits<TenT>::order_t;
    using BondLabelT = typename tensor_traits<TenT>::bond_label_t;
    OrderT num_bonds = a.Rank();
    OrderT num_cols = num_bonds - num_of_bonds_as_rows;
    List<BondLabelT> trans_labs(num_bonds);
    for(size_t k=0; k < static_cast<size_t>(num_cols); k++) {
      trans_labs[k] = num_of_bonds_as_rows + static_cast<OrderT>(k);
    }
    for(size_t k=static_cast<size_t>(num_cols); k < static_cast<size_t>(num_bonds); k++) {
      trans_labs[k] = static_cast<OrderT>(k) - num_cols;
    }
    auto adag = a;
    adag.Transpose(trans_labs);
    gqten::QR(&adag,num_cols,&q,&l);
    List<BondLabelT> q_labs(static_cast<size_t>(num_cols)+1);
    List<BondLabelT> l_labs(static_cast<size_t>(num_of_bonds_as_rows)+1);
    std::iota(q_labs.begin(),q_labs.end(),-1);
    q_labs[0] = num_cols;
    std::iota(l_labs.begin(),l_labs.end(),1);
    l_labs[static_cast<size_t>(num_of_bonds_as_rows)] = 0;
    l.Transpose(l_labs);
    q.Transpose(q_labs);
  }

  //----------------------------------------------------------------------------
  // eigh
  //----------------------------------------------------------------------------
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void eigh(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const order_t<TenT> &num_of_bonds_as_rows,
       real_ten_t<TenT> &w_diag,
       TenT &v) {
    using RealT = typename tensor_traits<TenT>::real_t;
    using ShapeT = typename tensor_traits<TenT>::shape_t;

    const auto shape_a = a.Shape();
    const size_t order_a = shape_a.size();
    const size_t ldims  = static_cast<size_t>(num_of_bonds_as_rows);
    assert(ldims > 0 && ldims <= order_a);

    size_t m = 1, ncols = 1;
    for (size_t i = 0; i < ldims; ++i)       m     *= shape_a[i];
    for (size_t i = ldims; i < order_a; ++i)  ncols *= shape_a[i];
    assert(m == ncols && "eigh: Hermitian matrix must be square after flattening");

    RealT * pw = nullptr;
    elem_t<TenT> * pv = nullptr;
    const char jobz = 'V';
    const char uplo = 'U';
    size_t n = 0;
    gqten::EigHerm(&a,ldims,pw,pv,&n,jobz,uplo);
    w_diag = gqten::tensor<RealT>({n},pw);
    ShapeT shape_v(static_cast<size_t>(num_of_bonds_as_rows)+1);
    for(size_t k=0; k < ldims; k++) {
      shape_v[k] = shape_a[k];
    }
    shape_v[ldims] = n;
    v = TenT(shape_v,pv);
  }

} // namespace tci

#endif // TCI_GQTEN_LINALG_H
