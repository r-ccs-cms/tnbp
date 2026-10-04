// main.cc
//
// Unit test for tnbp::SquareRootAndInverse (Hermitian eigh path) against a
// closed-form oracle, and for its agreement with the SVD-based reference
// implementation tnbp::SquareRootAndInverseViaSvd.
//
// Fixtures are M = W diag(lam) W^dagger with a deterministic unitary W
// (modified Gram-Schmidt of a fixed complex matrix; unitarity is asserted
// inside the test) and a known spectrum, so the expected R = W f(lam) W^dagger
// and S = W g(lam) W^dagger are computable in plain arrays independently of
// the code under test. Both the complex and the real tensor instantiations
// are exercised (the latter with a real orthogonal W).

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <type_traits>
#include <string>
#include <vector>

#include "tnbp/framework/typedef.h"
#include "tnbp/framework/root.h"

#include "typedef.h"

using RealTensor = typename tcapi::tensor_traits<Tensor>::real_ten_t;

using C = std::complex<double>;
using Mat = std::vector<std::vector<C>>;

static constexpr int N = 6;
static constexpr double TOL = 1.0e-10;

static int g_failures = 0;

void check(bool ok, const std::string & label, double value) {
  std::cout << (ok ? " PASS " : " FAIL ") << label
	    << " (value = " << value << ")" << std::endl;
  if (!ok) { g_failures++; }
}

Mat zeros_mat() { return Mat(N, std::vector<C>(N, C(0.0, 0.0))); }

Mat matmul(const Mat & a, const Mat & b) {
  Mat c = zeros_mat();
  for (int i = 0; i < N; i++)
    for (int k = 0; k < N; k++)
      for (int j = 0; j < N; j++)
	c[i][j] += a[i][k] * b[k][j];
  return c;
}

Mat adjoint(const Mat & a) {
  Mat c = zeros_mat();
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      c[i][j] = std::conj(a[j][i]);
  return c;
}

double max_abs_diff(const Mat & a, const Mat & b) {
  double d = 0.0;
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      d = std::max(d, std::abs(a[i][j] - b[i][j]));
  return d;
}

// Deterministic unitary (complex) or orthogonal (real) matrix via modified
// Gram-Schmidt of a fixed full-rank matrix. The entries grow as k^{3/2},
// which is nonlinear in the index and keeps every column linearly
// independent (an affine formula like i*N+j would give a rank-2 base whose
// later Gram-Schmidt columns are pure rounding noise). The columns are still
// nearly dependent, so a single Gram-Schmidt pass cannot reach
// machine-precision orthogonality; a second pass restores it, and main()
// asserts the result before any fixture relies on it.
Mat deterministic_unitary(bool complex_entries) {
  Mat a = zeros_mat();
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      a[i][j] = C(std::pow(double(i * N + j + 1), 1.5),
		  complex_entries ? std::pow(double(j * N + i + 1), 1.5)
				  : 0.0);
  for (int pass = 0; pass < 2; pass++) {
    for (int j = 0; j < N; j++) {
      for (int k = 0; k < j; k++) {
	C proj(0.0, 0.0);
	for (int i = 0; i < N; i++) proj += std::conj(a[i][k]) * a[i][j];
	for (int i = 0; i < N; i++) a[i][j] -= proj * a[i][k];
      }
      double nrm = 0.0;
      for (int i = 0; i < N; i++) nrm += std::norm(a[i][j]);
      nrm = std::sqrt(nrm);
      for (int i = 0; i < N; i++) a[i][j] /= nrm;
    }
  }
  return a;
}

// M = W diag(lam) W^dagger mapped through func on the spectrum.
Mat spectral_build(const Mat & w, const std::vector<double> & lam,
		   double (*func)(double)) {
  Mat c = zeros_mat();
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      for (int k = 0; k < N; k++)
	c[i][j] += w[i][k] * C(func(lam[k]), 0.0) * std::conj(w[j][k]);
  return c;
}

template <typename TenT>
typename tcapi::tensor_traits<TenT>::elem_t to_elem(const C & x) {
  using ElemT = typename tcapi::tensor_traits<TenT>::elem_t;
  if constexpr (std::is_floating_point_v<ElemT>) {
    return ElemT(x.real());
  } else {
    return ElemT(x);
  }
}

template <typename TenT>
TenT to_tensor(typename tcapi::tensor_traits<TenT>::context_handle_t & ctx,
	       const Mat & a) {
  using ShapeT = typename tcapi::tensor_traits<TenT>::shape_t;
  using CoorsT = typename tcapi::tensor_traits<TenT>::elem_coors_t;
  ShapeT shape(2);
  shape[0] = N;
  shape[1] = N;
  TenT t = tcapi::zeros<TenT>(ctx, shape);
  CoorsT coors(2);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      coors[0] = i;
      coors[1] = j;
      tcapi::set_elem(ctx, t, coors, to_elem<TenT>(a[i][j]));
    }
  }
  return t;
}

template <typename TenT>
Mat from_tensor(typename tcapi::tensor_traits<TenT>::context_handle_t & ctx,
		const TenT & t) {
  using CoorsT = typename tcapi::tensor_traits<TenT>::elem_coors_t;
  Mat a = zeros_mat();
  CoorsT coors(2);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      coors[0] = i;
      coors[1] = j;
      a[i][j] = C(tcapi::get_elem(ctx, t, coors));
    }
  }
  return a;
}

bool all_finite(const Mat & a) {
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      if (!std::isfinite(a[i][j].real()) || !std::isfinite(a[i][j].imag()))
	return false;
  return true;
}

template <typename TenT>
void run_fixture(typename tcapi::tensor_traits<TenT>::context_handle_t & ctx,
		 const std::string & name,
		 const Mat & w, const std::vector<double> & lam,
		 double sv_min, bool compare_with_svd_path) {
  std::cout << "=== fixture " << name << " ===" << std::endl;

  Mat m = spectral_build(w, lam, [](double x) { return x; });
  // Oracle threshold mirrors the documented contract: eigenvalues at or
  // below max(sv_min, 0) are dropped.
  double floor_ev = std::max(sv_min, 0.0);
  std::vector<double> lam_thr(lam);
  for (auto & x : lam_thr) { if (!(x > floor_ev)) x = 0.0; }

  Mat m_thr = spectral_build(w, lam_thr, [](double x) { return x; });
  Mat p = zeros_mat();
  for (int k = 0; k < N; k++)
    if (lam_thr[k] > 0.0)
      for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
	  p[i][j] += w[i][k] * std::conj(w[j][k]);

  TenT M = to_tensor<TenT>(ctx, m);
  TenT R;
  TenT S;
  tnbp::SquareRootAndInverse(ctx, M, R, S, sv_min);
  Mat r = from_tensor<TenT>(ctx, R);
  Mat s = from_tensor<TenT>(ctx, S);

  check(all_finite(r), name + ": R all finite", 0.0);
  check(all_finite(s), name + ": S all finite", 0.0);

  double herm = max_abs_diff(r, adjoint(r));
  check(herm <= TOL, name + ": R Hermitian", herm);

  double d_rr = max_abs_diff(matmul(r, r), m_thr);
  check(d_rr <= TOL, name + ": ||R.R - M_thr||", d_rr);

  double d_rs = max_abs_diff(matmul(r, s), p);
  check(d_rs <= TOL, name + ": ||R.S - P||", d_rs);

  double d_sms = max_abs_diff(matmul(s, matmul(m, s)), p);
  check(d_sms <= TOL, name + ": ||S.M.S - P||", d_sms);

  if (compare_with_svd_path) {
    TenT R_svd;
    TenT S_svd;
    tnbp::SquareRootAndInverseViaSvd(ctx, M, R_svd, S_svd, sv_min);
    double d_r = max_abs_diff(r, from_tensor<TenT>(ctx, R_svd));
    double d_s = max_abs_diff(s, from_tensor<TenT>(ctx, S_svd));
    check(d_r <= TOL, name + ": ||R_eigh - R_svd||", d_r);
    check(d_s <= TOL, name + ": ||S_eigh - S_svd||", d_s);
  }
}

double check_unitary(const Mat & w) {
  Mat eye_mat = zeros_mat();
  for (int i = 0; i < N; i++) eye_mat[i][i] = C(1.0, 0.0);
  return max_abs_diff(matmul(adjoint(w), w), eye_mat);
}

int main() {
  ContextHandle ctx;
  tcapi::create_context(ctx);
  typename tcapi::tensor_traits<RealTensor>::context_handle_t ctx_r;
  tcapi::create_context(ctx_r);

  // Self-verify the fixtures: W must be unitary / orthogonal to machine
  // precision, otherwise the spectral oracle is meaningless.
  Mat w = deterministic_unitary(true);
  double unit_dev = check_unitary(w);
  check(unit_dev <= 1.0e-12, "fixture: W unitary", unit_dev);
  Mat w_real = deterministic_unitary(false);
  double orth_dev = check_unitary(w_real);
  check(orth_dev <= 1.0e-12, "fixture: W_real orthogonal", orth_dev);

  // Fixture A: well-conditioned, every eigenvalue above the floor, so the
  // eigh path and the SVD reference path must agree.
  run_fixture<Tensor>(ctx, "A(well-conditioned)", w,
		      {1.0, 0.5, 0.1, 0.05, 0.01, 0.005}, 1.0e-8, true);

  // Fixture B: 17 decades of dynamic range with degenerate pairs and a tail
  // below the floor — the matrix class on which dense SVD misbehaves.
  run_fixture<Tensor>(ctx, "B(ill-conditioned)", w,
		      {1.0, 1.0, 1.0e-2, 1.0e-9, 1.0e-9, 1.0e-17}, 1.0e-8,
		      false);

  // Fixture C: a rounding-scale negative eigenvalue, the PSD-up-to-rounding
  // case the signed threshold must clamp to zero without producing NaN.
  run_fixture<Tensor>(ctx, "C(negative-tail)", w,
		      {1.0, 1.0e-2, 1.0e-9, -1.0e-13, 0.0, 0.0}, 1.0e-8,
		      false);

  // Fixture D: the real-tensor instantiation (real orthogonal W), covering
  // the RealTenT branch of the diagonal construction.
  run_fixture<RealTensor>(ctx_r, "D(real)", w_real,
			  {1.0, 0.5, 1.0e-2, 1.0e-9, 1.0e-13, 0.0}, 1.0e-8,
			  false);

  // Fixture E: negative sv_min, exercising the floor clamp max(sv_min, 0) —
  // every strictly positive eigenvalue is retained, none goes through sqrt
  // of a negative threshold.
  run_fixture<Tensor>(ctx, "E(negative-sv_min)", w,
		      {1.0, 0.5, 0.1, 0.05, 0.01, 0.005}, -1.0, false);

  std::cout << (g_failures == 0 ? " ALL PASS" : " FAILURES: ")
	    << (g_failures == 0 ? std::string("")
				: std::to_string(g_failures))
	    << std::endl;
  return g_failures == 0 ? 0 : 1;
}
