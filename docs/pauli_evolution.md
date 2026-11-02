# Pauli imaginary-time operators

The operator API constructs exponentials of a Hermitian Pauli-string sum on a
fixed simple graph. It does not run BP, change a state, measure energy, or print
logs. Include `pauli/evolution_plan.h` for tensor-independent planning and
`tnbp/operator/pauli_exponential.h` for TCAPI tensor construction. The latter is
also included by `tnbp/tnbp.h`.

1. `NormalizePauliHamiltonian(input)` accepts coefficient/string terms with
   complex-double coefficients. Labels use the existing convention: the
   rightmost character is site 0. Equal labels are merged; their resulting
   coefficients must be exactly real and finite. Identity terms are separated,
   exact zero terms are removed, and remaining terms are ordered by label.
2. `MakePauliEvolutionPlan(hamiltonian, ordered_edges)` includes every site
   `0..num_sites-1`, including isolated sites. New graph inputs use ascending
   endpoints while preserving edge order
   (see [graph conventions](graph_tensor_conventions.md)). The planner retains
   the supplied graph metadata for compatibility with existing states and
   tensor virtual axes, connects each term's support by a deterministic
   shortest-path tree, then first-fit packs pairwise commuting terms with
   disjoint tree sites into layers. Intermediate identity sites count as used
   sites, so even edge-disjoint terms split if their routes share a vertex. Disconnected multi-site support is rejected.
3. `MakeTrotterActions(plan, order=2)` returns one complete step in actual action
   order. Second order is the symmetric half-step/full-step/reverse-half-step
   sequence; zero layers yield no actions, one layer yields one full action.
4. `MakePauliTermExponentialTPO<T>(ctx, plan, term_id, delta)` and
   `MakePauliLayerExponentialTPO<T>(ctx, plan, layer_id, delta)` return
   `ScaledPauliTPO<T>{tensors, log_scale}`. These overloads cover all sites in site order.
   For each tensor the axes are incident edges in original edge-index order,
   then physical input and output. The value is `O[...,in,out] = <out|G|in>`.

Both generation APIs also accept a final `const std::vector<int>& target_sites`
argument. Only those tensors are constructed, in the supplied order; the list
may be empty or unsorted, but must contain unique valid sites. No ownership map
or MPI call is needed. `log_scale` is still for the complete term/layer and is
identical on all ranks; never sum it across ranks. An empty list returns no
tensors and the same whole-operator scale.

```cpp
auto local = tnbp::MakePauliLayerExponentialTPO<Tensor>(
    ctx, plan, layer_id, delta, site_idx);
tnbp::AttachTPO(ctx, plan.graph.edges, local.tensors,
    V, site_idx, owners, E, edge_idx, comm);
```

`site_idx` must match the local storage order of `V`. Plans, route metadata and
coefficients may be replicated, while tensor payloads and the CLI's layer cache
contain only local sites. The original all-site overloads remain available.

Plans are validated value objects and must not be mutated after construction.
Plans use version 2 for site-disjoint layering. Their versioned SHA256 digest covers the normalized problem and route/layer
plan, but excludes time steps and MPI ownership. No MPI or tensor backend is
needed for planning. Generation supports complex float/double tensor backends,
requires a finite nonnegative time step, and makes an identity TPO at zero time.
Each virtual dimension is at most two. A branching site's tensor size still
grows with its number of active virtual axes.

For a term, `x = delta * coefficient`, the exact exponential is
`exp(log_scale) * represented_TPO`, where `log_scale = abs(x)` and the represented
operator is `exp(-abs(x)) * exp(-x P)`. Stable bounded coefficients avoid directly
evaluating large cosh/sinh values. Layer scales add. Overflow of `x` or the
accumulated scale is rejected.

For the represented operator `a I + b P`, coefficients are shared equally
among the `m` non-identity Pauli sites: each contributes `a^(1/m)` on virtual
branch 0 and `abs(b)^(1/m)` on branch 1. A negative `b` contributes its sign
only at the deterministic Pauli root. Each site is weighted once, including
Pauli sites with several tree legs. Intermediate I sites and I branch points
remain unweighted copy tensors. The virtual flags still agree along the tree,
so their products are `a` and `b` up to rounding. Single-site terms retain
their direct local exponential; the common `log_scale` is unchanged.

This layout is identified by `PauliExponentialWeights` (`pauli-support-equal-v1`).
It balances the virtual representation without changing the exact operator or
the root/inverse cutoff and compression algorithms. Finite cutoff projections
are representation-dependent; balancing does not guarantee globally exact
compression or monotonic energy on a loopy graph.

The separated identity Hamiltonian coefficient is **not** included in any term
or layer TPO. The caller accounts for its log factor `-time_step * c_identity`
exactly once per complete Trotter step. Operator scales are not the norm of a
subsequently compressed state.

All-site tensors can be consumed by `AbsorbTPO`: the explicit-site overload
accepts isolated sites, while the legacy overload infers sites from edges.
Local tensors use `AttachTPO`, including isolated sites. The CLI normalizes
vertices not visited by compression locally after each action. See
[tests](../test/pauli_evolution/README.md) for exact references and validation
scope.

The imaginary-time CLI is in `apps/imaginary_time_evolution`. Its existing
BP-based compression is approximate: local SVD tails do not include losses from
inverse-square-root support projections. Operator-construction exactness and
Trotter order checks are separate from that compression approximation.

For caller-controlled BP between groups of independent edges, see the optional
[batch compression API](truncation_batches.md) and
`MakePauliLayerTruncationBatches(plan, layer_id)`. The CLI uses automatic
leaf-side batches with BP between them; the all-edge API remains available.
