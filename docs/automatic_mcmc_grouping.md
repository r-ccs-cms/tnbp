# Automatic three-line grouping

`make_mcmc_grouping(edges, lines)` in `mcmc_plan.h` constructs the same
`MCMCGrouping` accepted by manual input. It reads the original graph without
contracting it or reordering its edges. For each target line, it partitions the
sites in the immediately lower and upper lines into groups anchored at the
ordered target sites. End lines have only one neighboring line.

## Selection policy

Each neighboring line is partitioned independently:

1. Visit target anchors in their specified line order. Reserve one as-yet
   unassigned site connected by an original edge, if available. Ties use the
   neighboring line's specified site order. A prior anchor cannot consume all
   direct neighbors before later anchors get this first opportunity.
2. Grow the groups over original edges. Available candidates are direct
   neighbors of an anchor or neighbors of sites already assigned to that group.
   Process the shortest path length from the anchor first, then anchor position,
   then candidate position in the neighboring line. Never select an already
   assigned site or reassign it to another group.
3. Reject an unassigned component with no path to an anchor. Singleton target
   groups are allowed. All neighboring sites must be covered exactly once.

Distance counts unweighted original edges along these admissible paths. A path
may use the anchor and sites in the neighboring line that belong to its group;
it cannot pass through another target anchor or a site owned by another group.
This constraint makes each selected member absorbable after its predecessors.
It is not a shortest path through arbitrary sites outside the window. The initial
reservation has priority over later growth; this is a deterministic greedy
partition, not a global matching or partition optimizer.

Each final group lists its target first, then lower-line members in assignment
order, then upper-line members in assignment order. Ties depend on line order,
not numeric site labels, edge-vector order or MPI rank. Reversing edge endpoints
also leaves the grouping unchanged. Tensor axes still follow edge-vector order;
see [graph conventions](graph_tensor_conventions.md).

For two chains with rungs 0--6, 2--8 and 4--10:

```text
lines[0] = {0,1,2,3,4,5}
lines[1] = {6,7,8,9,10,11}
groups[1] = {{6,0,1}, {7}, {8,2,3}, {9}, {10,4,5}, {11}}
```

The non-rung sites are covered by growth, rather than silently omitted.

## One plan for all execution paths

Initialization resolves automatic input to a concrete plan before constructing
boundaries. The workspace builds it once across batches; `site_line_sets` stores
the resolved plan. Sampling reuses that plan for both sweep directions. It does
not reconstruct groups from compressed boundary graphs.

Automatic and manual plans use the same validation, graph-only increments,
tensor absorption, post-absorption bond relocation, and left/right environment
construction. The standalone two-line BMPS API without explicit groups uses the
same neighbor-partition primitive; callers implementing repeated sweeps should
pass a fixed plan, as MCMC does.

`manual_grouping` is an input override: empty means automatic initialization.
Changing a nonempty override after initialization is rejected if it disagrees
with `site_line_sets`. To construct a plan, use `make_mcmc_grouping` with the
original graph and all lines. A contracted graph and one line are insufficient
to reconstruct the original three-line plan.

## Supported layouts and validation

Both modes require at least two nonempty lines covering the graph exactly once,
with edges only within the same or adjacent lines. Each group and each left/right
sweep prefix/suffix must be connected. Planned absorptions must remain possible
on the evolving graph, including preceding bond relocations.

The automatically generated candidate is checked with `validate_mcmc_grouping`
before numerical execution. If it fails, the call raises `invalid_argument`;
it never silently switches groups during a sweep. A different manual partition
or line layout may work. This policy does not guarantee a valid partition for
every graph and every line order.
