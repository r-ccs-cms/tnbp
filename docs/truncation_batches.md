# Batched bond compression

`tnbp::TruncationBatch` compresses a caller-selected set of edges using the
current BP messages. It does not run BP. The existing `tnbp::Truncation` API
still processes the complete graph in its original edge order. The
imaginary-time CLI now uses the leaf-side schedule described below.

Include `tnbp/ptns/function/truncation.h` or `tnbp/tnbp.h`.

## Selection and results

Pass global edge IDs: positions in the complete ordered `edges` vector, not
positions in a rank's `edge_idx`. The selected edges must share no sites.
Unordered input is accepted and processed in increasing edge-ID order.
Duplicates, out-of-range IDs and edges sharing an endpoint are rejected before
numerical updates. Participating ranks also reject inconsistent local tensor/index
counts. These checks are local: callers must provide compatible valid inputs on
both owners, since a failure on just one owner cannot cancel its peer. Empty
selection is a no-op.

```cpp
auto result = tnbp::TruncationBatch(
    ctx, edges, batch_edge_ids,
    V, site_idx, owners, E, edge_idx, comm,
    max_bond_dim, sv_min, truncation_error);
```

`TruncationBatchResult<Tensor>` contains parallel arrays:

- `edge_ids`: selected edges incident to this rank, in global edge-ID order.
- `bond_dims`: retained dimensions for those edges.
- `truncation_errors`: their existing per-edge SVD error estimates.

There are no entries for unselected edges. A boundary edge returns a result on
both endpoint owners; count it only on the owner of `edges[e].first` when forming
global statistics. Optional final argument `site_norms` retains the existing
append-only normalization accounting for selected locally incident edges.

Only endpoint tensors and both directed messages of selected edges are changed.
Other tensors/messages retain their values and bytes. Tensor axes, graph order,
site ownership and local index arrays are unchanged. Numerical decomposition,
cutoffs and the authoritative first-endpoint-owner protocol are shared with the
all-edge API.

## Two independent three-site operations

For `0—1—2` and `3—4—5`, let the ordered graph be
`{{0,1},{1,2},{3,4},{4,5}}`. The first bonds have IDs `{0,2}` and the second bonds
have IDs `{1,3}`. Both are valid batches:

```cpp
auto first = tnbp::TruncationBatch(
    ctx, edges, {0,2}, V, site_idx, owners, E, edge_idx, comm,
    max_bond_dim, sv_min, truncation_error);

// Optional: run the caller's existing BP convergence loop here, using
// UpdateBPMessages / ExchangeBPMessages / BeliefPropagationDifference.

auto second = tnbp::TruncationBatch(
    ctx, edges, {1,3}, V, site_idx, owners, E, edge_idx, comm,
    max_bond_dim, sv_min, truncation_error);
```

The optional BP loop belongs to the caller. The batch API neither adds an
implicit convergence step nor requires a particular batching strategy. A batch
can also contain a single edge. Site disjointness avoids conflicting updates;
it is not a guarantee of exact compression on an arbitrary state or graph.

## MPI contract

The complete graph, owners, selection and numerical settings must already be
agreed by the caller. Validation is local; the API adds no collective agreement
or barrier. Both owners of a selected edge must enter with the same batch.
Other ranks can call and return an empty result if uninvolved, or omit the call.
As with the all-edge routine, backend/communication failures after calculation
starts may leave partial updates.

Independent rank pairs can progress concurrently. Edges owned by the same rank
are processed locally in sequence. No new tensor gather, communicator, OpenMP
threading or global synchronization is introduced. This preserves the existing
owner-local communication model while allowing BP between batches when desired.

See [batch tests](../test/pauli_evolution/README.md#batch-compression-api) for
six-site exact references and MPI 1/2/3/4-rank checks.

## Automatic leaf-side schedules for Pauli layers

`MakePauliLayerTruncationBatches(plan, layer_id)` returns ordered batches of
global edge IDs. It does not run MPI, modify the graph, or compress tensors.

Each term's route is rooted at its existing `root_site`. An edge becomes ready
only after all edges below its child endpoint have been processed. The smallest
global edge ID among ready edges is taken next. This gives a deterministic
leaf-to-root sequence; it does not collect multiple edges of one tree into a
batch, even when those edges share no endpoint.

Visit trees in increasing term ID order. Only each tree's next edge is eligible.
Append it if neither endpoint is occupied in the batch; otherwise leave that
tree waiting without skipping ahead. Repeat until all sequences are exhausted.
Every batch contains at most one edge per tree and has no shared endpoints.
Single-site terms have empty sequences. MPI owners do not affect this schedule.
The default planner now puts site-disjoint trees in a layer, including all
intermediate identity sites. Its tree heads therefore do not conflict; the
batch builder retains the endpoint guard as a defensive check.

For a branching tree with `{{0,1},{0,2},{1,3},{1,4}}`, the sequence is
`{1,2,3,0}`. For two trees with edges
`{{0,1},{1,2},{1,3},{4,5},{5,6}}` and respective support roots 0 and 4, the
sequences are `{1,2,0}` and `{4,3}`, giving batches `{1,4}`, `{2,3}`, `{0}`.

The CLI performs `absorb -> BP -> batch 0 -> BP -> batch 1 -> ...`.
Only layer-used edges are compressed. There is no BP after the last batch and
no extra measurement BP. BP limits/tolerances remain the configured values.
Vertices not visited by compression, including connected single-site terms,
are normalized locally. The all-edge API and batch numerical kernel are
unchanged. Leaf order is a fixed choice, not a general exactness guarantee.
