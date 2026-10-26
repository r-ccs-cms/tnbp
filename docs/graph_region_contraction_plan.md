# Small region contraction plans

Include `tnbp/framework/graph_region_contraction_plan.h`. These free functions
need no tensor backend or MPI. A region is a sorted list of unfused edge tokens;
each token has a positive dimension and occurs in at most two input regions.
A token appearing once is an external leg. Outer products are not considered.

`make_small_region_contraction_plan(leaf_axes, token_dimensions, leaf_owners,
repeated_leaf, repetitions, fallback_steps, steps, max_elements, used_fallback)`
returns total dense scalar MACs. For up to five leaves it searches every
connected binary tree (at most 180 pair-choice sequences). Above five it checks
and uses the supplied fallback. The fallback must itself be a connected tree.

A step is `{lhs, rhs, result}` with leaf IDs `0..n-1` and result ID `n+k` for
step `k`. This is **not** the edge schedule accepted by
`graph_tensor_network_contraction`. Every non-root node is consumed once.
The operand containing `repeated_leaf` stays on the left; otherwise the operand
containing the smallest leaf ID stays on the left. The result inherits that
owner. For more than five leaves the fallback's operand directions are retained,
while the repeated leaf must still remain on the left.

The objective counts common subtrees once and repeated-dependent steps
`repetitions` times. Ties use largest output element count, then normalized
steps in lexicographic order. Common steps precede dependent steps. This does
not optimize communication, transpose/copy traffic, wall time or peak live
memory. The maximum output size is zero for a singleton and one for a scalar.

`describe_region_contraction_plan` derives node axes, owners, dependency flags,
parents, common frontier and per-step costs/sizes. The frontier contains the
common nodes borrowed by repeated steps; other temporary nodes can be released
after their parent executes. A plan describes shapes, not reusable tensor values.

`make_region_join_metadata` creates the same two-operand graph and fusion
instructions for both cost evaluation and numerical execution. Shared tokens
are fused in ascending order; each external leg gets a separate metadata leaf.
The cost is evaluated with `graph_tensor_contraction_cost`. Numerical execution
uses the existing distributed `graph_tensor_contraction`; external metadata
leaves require no tensors or owners.

Invalid input throws. An overflowing candidate is excluded from the small
search; failure of every candidate or overflow in a large fallback throws.
Resource failures are not treated as rejected candidates. All output arguments
must be distinct and remain unchanged on failure.

These token-based functions remain standalone helpers and regression references.
The default full-state heat-bath path now uses original-edge symbolic plans;
see [Symbolic graph contraction costs](graph_symbolic_contraction_cost.md).
Its numerical executor and communicator lifetime are described there. Legacy
sampling methods are unchanged.
