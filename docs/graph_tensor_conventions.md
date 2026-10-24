# Graph, tensor-axis and ownership conventions

These conventions apply to the operations in
[`graph_contraction.h`](../include/tnbp/framework/graph_contraction.h), including
the graph bond-shift functions used by boundary MPS.

## Edges determine tensor-axis order

`edges` is an ordered vector of undirected site pairs. For a site, scan this
vector from beginning to end and keep its incident edges. Tensor axis `k`
corresponds to the `k`-th edge in that filtered list, as returned by
`GetSurroundingBondIndex`. Site-label order does not determine axis order.

For example:

```cpp
edges = {{2, 5}, {0, 2}, {2, 7}, {0, 5}};
// Tensor at site 2: axes correspond to neighbors [5, 0, 7].
// Tensor at site 5: axes correspond to neighbors [2, 0].
```

The dimension of a bond must agree at both tensor endpoints. Open or physical
legs in these graph operations must also be represented by edges; do not append
implicit tensor axes that have no corresponding edge.

Reversing a pair from `(a,b)` to `(b,a)` does not change its meaning or axis
position. Reordering the vector itself does change axis order and requires
matching tensor transpositions. Do not sort `edges` independently of tensors.

Graph operations may replace, remove or append edges. They update affected
tensor axes along with the graph, so **on return, axis order follows the updated
`edges` vector**. They do not promise to preserve old edge indices. Recompute
cached incident-edge indices after a graph update.

For bond shifts, retained source axes keep their relative edge order and the
new source--destination SVD bond is appended. Destination and shared-neighbor
axes are transposed/reshaped as needed by graph contraction and bond fusion.

## Site storage and owners

The ownership map has the meaning:

```cpp
std::map<IntT, IntT> owners; // site label -> MPI rank in comm
// On the rank owners.at(site_idx[k]):
// tensors[k] is the tensor at site_idx[k].
```

Map keys are site labels, not positions in `tensors`, `site_idx` or `edges`.
The map's sorted iteration order defines neither tensor storage nor tensor-axis
order. Local storage is determined only by the parallel vectors `site_idx` and
`tensors`; every supplied local site must belong to the calling rank.
Rank values refer to the supplied communicator `comm`.

Site labels need not be consecutive or start at zero. Generic graph operations
and portable checkpoints use signed integer labels; a label is not a dense
array index. Portable/app interfaces here use `int` labels. Rank values are
separate and must be nonnegative communicator ranks. The QASM apps use `q[i]`
indices as labels, so their circuit operands must be nonnegative even though
the graph/checkpoint representation can hold negative labels.


All participating ranks must agree on graph order, operation arguments and
ownership metadata. Tensor data are stored on their owners. A graph can include
sites whose tensors are absent from the supplied vectors, provided the operation
does not need to contract or reshape those tensors. Required owner entries must
be present; missing ownership is not a default to rank zero.

The main estimator and sampling apps accept an optional ownership file and
print the selected map with `#`-prefixed lines. Initialization consumes that
map directly; it does not choose a partition internally. See
[application site ownership](site_ownership.md).

## Ownership for each operation

| Operation | Result ownership | Required tensors |
| --- | --- | --- |
| `graph_tensor_contraction(..., site_a, site_b, site_c, comm)` | `site_c` inherits the original owner of `site_a`; `site_b` is removed | `site_a`, `site_b`, and their common neighbors, whose bonds are fused |
| `graph_shift_bond_to_neighbor(..., source, destination, other, comm)` | All original sites retain their original owners | `source`, `destination`, and `other` only if destination--other already exists |
| `graph_shift_bonds_to_neighbor(..., source, destination, others, comm)` | All original sites retain their original owners | `source`, `destination`, and selected others with an existing destination--other edge |

For contraction, `site_c` must be `site_a` or a fresh label. The current function
does not support `site_c == site_b`. Swapping the operands changes the owner of
the result. Unaffected graph sites need not have tensors or owner entries;
every supplied local tensor still needs an owner entry.

## Bond-shift contract

```cpp
tnbp::graph_shift_bond_to_neighbor(
    ctx, edges, tensors, site_idx, owners, source, destination, other, comm);
```

This moves source--other to destination--other and always leaves a new
source--destination SVD bond. The three sites must be distinct and the selected
source--other edge must exist. The graph must be simple (no self-loops or
duplicate undirected edges).

- If source--destination already exists, it is included in the SVD factor and
  replaced by the new SVD bond.
- If it does not exist, a temporary dimension-one edge permits absorption as an
  outer product. Initial adjacency is therefore not required.
- If destination--other exists, it is fused with the moved bond at both ends.
- Unselected source bonds stay at source. Site labels, owners and local site
  order are preserved, while tensor values, shapes and edges are updated.

Internally, the temporary SVD factor belongs to the source owner and is absorbed
on the destination owner, with MPI transfer when necessary. Its temporary site
and ownership entry are removed. **Moving a bond does not transfer ownership of
either original site.**

The plural function takes a nonempty, unique vector of other endpoints and moves
their bonds together in one SVD. Both forms use full TCAPI SVD without truncation
or normalization; bond dimensions may grow. The generic functions do not choose
a line order. `bmps_shift_bonds_to_neighbors` chooses moves using the supplied
line order, then delegates to the plural function. It runs automatically after absorption
in `bmps_increment_line`, before subsequent BP/truncation.

`graph_shift_bonds_to_neighbor_for_edges` reproduces the ordered edge update
without tensors. `bmps_increment_line_for_edges` applies both absorption and
these shifts, using the same selection loop as the tensor path. MCMC boundary
plans, grouping validation and left/right environments use these updated graphs.
The current target line is not independently shifted in the environment path:
its neighboring boundary tensors already follow their post-absorption graphs.

Structural errors are rejected before mutation under the shared-input contract.
A backend failure during final absorption can leave partial output, as with
`graph_tensor_contraction`; these operations are not transactional.

## Validation

See [bond-shift tests](../test/graph_tensor_contraction/README.md). CPU tests cover real and
complex tensors, MPI 1/2/3 ranks, reversed edge order, partial tensor/owner
inputs, shared-neighbor fusion and absent neighbor bonds. They compare every
open-boundary tensor component before/after and check site/owner preservation.
CUDA tests have passed on one GPU and same-node two GPUs, with one visible GPU
per MPI process. Three-GPU execution has not been validated.

## Shared validation and bond fusion

`detail::fuse_graph_bond_axes` implements the common transpose/reshape for the
contracted tensor and each shared neighbor. Its arguments specify the kept and
removed axes: the fused dimension has the kept axis followed by the removed
axis, while every other axis keeps its relative order. This works whether the
removed edge originally precedes or follows the kept edge.

The public `graph_tensor_contraction` uses the same collective
`detail::validate_graph_tensors` preflight as bond shifts and BMPS sweeps.
It rejects nonadjacent/identical operands, result-label collisions (including
external graph labels or isolated owner entries), self/duplicate edges,
missing required tensors/owners, invalid ranks, inconsistent local storage,
wrong required tensor axis counts and incompatible internal bond dimensions.
Rejection happens on all ranks before point-to-point communication or mutation.
Graph order, operation arguments and owner metadata must still agree on ranks;
this preflight does not compare independent rank-local copies of that metadata.

Sampling also has an explicitly internal **partial environment** path,
`detail::graph_tensor_network_partial_contraction`. Here owner `-1` denotes a
graph-only site. Each operand pair must be either both present or both
graph-only. Present shared-neighbor tensors are fused normally, while absent
boundary tensors are reformed separately by the sampling caller before the
final environment contraction. This supports both boundary-axis reform without
operand tensors and local contraction with an absent complementary environment.
It is not an exception to the public contraction contract: public contraction
and bond shifting reject required owner `-1` entries. Supplied local tensors
always need valid owners, even in the partial path.

Preflight adds MPI collective operations to each contraction. It does not make
backend operations transactional or guarantee that a later allocation/SVD/
contraction failure leaves all tensors unchanged.

## One edge-update rule for contraction

`detail::make_graph_contraction_plan` constructs an immutable plan from the
original ordered graph. It supplies operand incident-edge indices, the
unfused result-axis order, ordered shared-neighbor fusion records, and the
updated edge vector. Original edge indices remain valid throughout tensor
fusion; the graph is replaced only after those operations.

Both tensor contraction (including partial-environment reform) and
`graph_tensor_contraction_for_edges` use this plan. The graph-only path performs
no TCAPI operation or MPI communication. Its rules are:

1. Remove the contracted A--B edge.
2. For a shared neighbor P, keep A--P in its original relative position and
   remove B--P; fuse axes in A-edge, B-edge order at both tensor endpoints.
3. Relabel B-only edges in place, rather than appending them.
4. Relabel surviving A/B endpoints to C and normalize those endpoint pairs.
5. Preserve the relative order and endpoint representation of unaffected edges.

The plan rejects invalid/simple-graph violations, absent or identical operands,
C==B and result-label collisions in the graph before returning. The tensor path
also validates owner-map collisions and tensor requirements collectively.
The graph-only function leaves its input unchanged when planning fails.

Bond shifting shares `detail::make_graph_shift_edges` for temporary-factor edge
construction, then uses this same contraction plan when absorbing the factor.
Temporary-label selection still accounts for owner-only sites in the tensor
path; different collision-free temporary labels do not survive absorption.
