// main.cc
//
// Unit test for tnbp::Truncation's retained bond dimension and reported
// truncation error against a closed-form oracle, on a single rank.
//
// The fixture is the smallest graph Truncation accepts — two sites joined by
// one edge — with both messenger tensors set to E = diag(lambda) for a
// strictly positive lambda, and sv_min = 0 so SquareRootAndInverse's
// max(sv_min, 0) floor drops nothing. Then Ra = Rb = diag(sqrt(lambda)) as a
// matrix function (basis-invariant, so the eigensolver's ordering and phases
// do not enter), and the contraction Truncation performs,
//
//   T[a,b] = sum_k Ra[k,a] Rb[k,b] = delta_ab * lambda_a,
//
// makes T = diag(lambda) up to the uniform normalization that follows. Its
// singular values are therefore lambda itself, known exactly rather than
// measured, so the retained chi and the reported error are both predictable
// in closed form from lambda alone.
//
// The rule predicted is the TCAPI v1 one: retain the smallest chi in
// [chi_min, chi_max] whose relative truncation error
//
//   epsilon(chi) = sum_{i >= chi} s_i^2 / sum_i s_i^2
//
// is at most the target. That is the spec's contract, asserted here in place
// of any particular backend's reading of the argument.
//
// The two readings of "at most" — strict and non-strict — diverge only where
// an epsilon equals the target exactly, so no case may sit there. Each case
// checks that for its own target rather than resting on a claim about lambda,
// which would go stale the moment either changes.
//
// Cases B and C differ in nothing but the target error, so a build in which
// Truncation ignores that argument fails C against B's outcome.

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "mpi.h"

#include "tnbp/tnbp.h"

#include "typedef.h"

static constexpr double TOL = 1.0e-12;

static const std::vector<double> kLambda =
  { 1.0, 0.5, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4 };

static int g_failures = 0;

void check(bool ok, const std::string & label, double value) {
  std::cout << (ok ? " PASS " : " FAIL ") << label
	    << " (value = " << value << ")" << std::endl;
  if (!ok) { g_failures++; }
}

// sum_{i >= chi} lambda_i^2, summed over exactly the values chi discards.
double tail_weight(int chi) {
  double t = 0.0;
  for(std::size_t i = static_cast<std::size_t>(chi); i < kLambda.size(); i++) {
    t += kLambda[i] * kLambda[i];
  }
  return t;
}

double epsilon_at(int chi) { return tail_weight(chi) / tail_weight(0); }

// The chi the TCAPI v1 rule retains. epsilon is decreasing in chi, so the
// first chi meeting the bound is the smallest one that does; chi_max then
// caps it and chi_min floors it.
int oracle_chi(int chi_min, int chi_max, double err) {
  const int k = static_cast<int>(kLambda.size());
  int chi_target = k;
  for(int chi = chi_min; chi <= k; chi++) {
    if( epsilon_at(chi) <= err ) { chi_target = chi; break; }
  }
  return std::max(chi_min, std::min(chi_target, chi_max));
}

// E = diag(lambda), the messenger tensor both ends of the edge carry.
Tensor make_messenger(ContextHandle & ctx) {
  using ShapeT = typename tcapi::tensor_traits<Tensor>::shape_t;
  using CoorsT = typename tcapi::tensor_traits<Tensor>::elem_coors_t;
  const int n = static_cast<int>(kLambda.size());
  ShapeT shape(2);
  shape[0] = n;
  shape[1] = n;
  Tensor e = tcapi::zeros<Tensor>(ctx, shape);
  CoorsT coors(2);
  for(int i = 0; i < n; i++) {
    coors[0] = i;
    coors[1] = i;
    tcapi::set_elem(ctx, e, coors, Elem(kLambda[i]));
  }
  return e;
}

// A site tensor with one leg per surrounding bond, per the convention
// GetSurroundingBondIndex fixes; the two-site graph gives each site exactly
// one. The entries only have to be nonzero, since Truncation renormalizes the
// site after applying its factor.
Tensor make_site(ContextHandle & ctx) {
  using ShapeT = typename tcapi::tensor_traits<Tensor>::shape_t;
  using CoorsT = typename tcapi::tensor_traits<Tensor>::elem_coors_t;
  const int n = static_cast<int>(kLambda.size());
  ShapeT shape(1);
  shape[0] = n;
  Tensor v = tcapi::zeros<Tensor>(ctx, shape);
  CoorsT coors(1);
  for(int i = 0; i < n; i++) {
    coors[0] = i;
    tcapi::set_elem(ctx, v, coors, Elem(double(i + 1)));
  }
  return v;
}

// Truncation rewrites V and E in place, so each case rebuilds both.
void run_case(ContextHandle & ctx, MPI_Comm comm,
	      const std::string & name, int max_dim, double err) {
  std::cout << "=== case " << name << " (max_dim = " << max_dim
	    << ", err = " << err << ") ===" << std::endl;

  const std::vector<std::pair<int,int>> I = { {0,1} };
  const std::vector<int> site_idx = { 0, 1 };
  const std::map<int,int> site_to_mpi_rank = { {0,0}, {1,0} };
  const std::vector<int> edge_idx = { 0 };

  // A target sitting on one of the epsilons is the one input on which the
  // strict and non-strict readings of the bound disagree, and the expectation
  // below would then encode whichever reading the backend happened to take.
  // A zero target is exempt: it says "no target truncation" to the spec and to
  // oracle_chi alike, rather than naming a bound some chi could sit on.
  if( err > 0.0 ) {
    double closest = 1.0;
    for(int chi = 0; chi <= static_cast<int>(kLambda.size()); chi++) {
      closest = std::min(closest, std::abs(epsilon_at(chi) / err - 1.0));
    }
    check(closest > 1.0e-6, name + ": target clear of every epsilon", closest);
  }

  std::vector<Tensor> V = { make_site(ctx), make_site(ctx) };
  std::vector<Tensor> E = { make_messenger(ctx), make_messenger(ctx) };

  std::vector<BondDim> res_bond_dim;
  std::vector<Real> res_trunc_err;

  tnbp::Truncation(ctx, I, V, site_idx, site_to_mpi_rank, E, edge_idx, comm,
		   static_cast<BondDim>(max_dim),
		   static_cast<Real>(0.0),
		   static_cast<Real>(err),
		   res_bond_dim, res_trunc_err);

  const int chi = oracle_chi(1, max_dim, err);
  const double expected_err = epsilon_at(chi);

  check(res_bond_dim.size() == 1 && res_trunc_err.size() == 1,
	name + ": one result per edge",
	double(res_bond_dim.size()));
  if( res_bond_dim.size() != 1 || res_trunc_err.size() != 1 ) { return; }

  check(static_cast<int>(res_bond_dim[0]) == chi,
	name + ": retained bond dim (expected " + std::to_string(chi) + ")",
	double(res_bond_dim[0]));

  const double d_err = std::abs(double(res_trunc_err[0]) - expected_err);
  check(d_err <= TOL * std::max(1.0, expected_err),
	name + ": reported trunc_err (expected "
	+ std::to_string(expected_err) + ")",
	double(res_trunc_err[0]));

  // The retained chi is what consumers actually feel: both incident sites and
  // both messenger slots must carry the truncated bond, not just the report.
  auto shape_v0 = tcapi::shape(ctx, V[0]);
  auto shape_v1 = tcapi::shape(ctx, V[1]);
  auto shape_e0 = tcapi::shape(ctx, E[0]);
  auto shape_e1 = tcapi::shape(ctx, E[1]);
  check(static_cast<int>(shape_v0[0]) == chi, name + ": V[0] bond", double(shape_v0[0]));
  check(static_cast<int>(shape_v1[0]) == chi, name + ": V[1] bond", double(shape_v1[0]));
  check(static_cast<int>(shape_e0[0]) == chi && static_cast<int>(shape_e0[1]) == chi,
	name + ": E[0] square of chi", double(shape_e0[0]));
  check(static_cast<int>(shape_e1[0]) == chi && static_cast<int>(shape_e1[1]) == chi,
	name + ": E[1] square of chi", double(shape_e1[0]));
}

int main(int argc, char * argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int mpi_size; MPI_Comm_size(comm, &mpi_size);
  if( mpi_size != 1 ) {
    std::cerr << " this test is single-rank; run with -np 1 " << std::endl;
    MPI_Abort(comm, 1);
  }

  ContextHandle ctx;
  tcapi::create_context(ctx);

  // Self-verify the fixture: the oracle below assumes lambda is strictly
  // positive and strictly decreasing, which is what makes epsilon monotone and
  // the SVD's descending order equal to lambda's own.
  bool spectrum_ok = kLambda.back() > 0.0;
  for(std::size_t i = 1; i < kLambda.size(); i++) {
    if( !(kLambda[i] < kLambda[i-1]) ) { spectrum_ok = false; }
  }
  check(spectrum_ok, "fixture: lambda positive and strictly decreasing",
	kLambda.back());

  // Case A: chi_max binds alone. With no target error this is the behaviour
  // the pre-restore call had, so A pins that the restore left it intact.
  run_case(ctx, comm, "A(chi_max-only)", 4, 0.0);

  // Case B: neither constraint binds — the full spectrum survives.
  run_case(ctx, comm, "B(no-truncation)", 6, 0.0);

  // Case C: B with a target error that admits chi = 2, which is the whole
  // capability the err argument names. epsilon(2) = 8.0e-3, epsilon(1) = 0.21.
  run_case(ctx, comm, "C(target-error)", 6, 1.0e-2);

  // Case D: the same target error under a chi_max that is tighter than what
  // the error alone would keep, so the cap wins and the reported error grows
  // to match what the cap actually discarded.
  run_case(ctx, comm, "D(chi_max-beats-target)", 1, 1.0e-2);

  // Case E: the bound alone drives chi down to one, with nothing else
  // reaching it — chi_max is slack and epsilon(1) = 0.21 <= 0.5 < epsilon(0).
  run_case(ctx, comm, "E(bound-drives-to-one)", 6, 0.5);

  // Case F: a target no epsilon can exceed, since epsilon(0) = 1. The bound
  // would admit discarding the whole spectrum, so the chi_min = 1 floor is the
  // only thing left deciding the answer — which E cannot show, its target
  // being one epsilon(0) already fails.
  run_case(ctx, comm, "F(chi_min-floor)", 6, 1.5);

  std::cout << (g_failures == 0 ? " ALL PASS" : " FAILURES: ")
	    << (g_failures == 0 ? std::string("")
				: std::to_string(g_failures))
	    << std::endl;

  MPI_Finalize();
  return g_failures == 0 ? 0 : 1;
}
