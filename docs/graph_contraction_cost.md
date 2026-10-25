# Graph contraction cost estimates

Include `tnbp/framework/graph_contraction_cost.h`. The evaluator uses ordered
edges and global site shapes; it needs no tensor backend, MPI initialization,
or tensor allocation. It evaluates a supplied schedule, without searching for
an order or changing MCMC behavior.

## Example

```cpp
std::vector<std::pair<int,int>> edges{{0,1},{0,2},{1,2}};
std::map<int,std::vector<std::uint64_t>> shapes{
    {0,{2,3}}, {1,{2,5}}, {2,{3,5}}};
auto schedule = edges;
std::vector<std::uint64_t> macs, elements;
auto total = tnbp::graph_tensor_network_contraction_cost(
    edges, shapes, schedule, {}, {}, &macs, &elements);
// total == 45; macs == {30,15}; elements == {15,1}
// edges is empty; shapes contains site 0 with scalar shape {}.
```

Each site's axes must follow its incident-edge order in `edges`. Supply shapes
for every graph vertex, including neighbors whose axes will be fused. This is
complete global metadata, rather than the numerical executor's rank-local tensor
list. Extra isolated scalar sites may be supplied with shape `{}`. Dimensions
must be positive integral values, and the dimensions at each bond's endpoints
must agree.

## Functions

- `graph_tensor_contraction_cost(edges, shapes, site_a, site_b, site_c)`
  evaluates one binary contraction and returns its dense MAC count.
- `graph_tensor_network_contraction_cost(edges, shapes, schedule,
  site_priority, label_map, ...)` evaluates the edge schedule with the same
  priority, tie-breaking, relabeling, and schedule updates as
  `graph_tensor_network_contraction` and returns total dense MACs. The last two
  input maps default to empty.

Both functions update their `edges` and `shapes` arguments. Pass independent
copies of the same initial metadata when comparing schedules. Partial and empty
schedules are supported. Shared graph plans update common-neighbor shapes as
well as the result shape; edge fusion and surviving edge order follow the
numerical API's rules. `site_c` must be `site_a` or a fresh label.

The network evaluator optionally fills these output containers, in this order:

| Output pointer | Per executed binary contraction |
|---|---|
| `step_macs` | Dense scalar multiply-accumulates |
| `step_elements` | Result tensor's element count |
| `steps` | `{site_a, site_b, site_c}` |
| `shape_history` | All remaining site shapes after the step |
| `edge_history` | Ordered edges after the step |

Pass `nullptr` for reports you do not need. Histories have no initial entry.
A schedule's redundant edges can disappear through fusion, so report length
can be smaller than the original schedule length. Optional output containers
must be distinct. Validation or integer overflow leaves inputs and reports
unchanged.

## Cost meaning

For one contraction, the dense MAC count is the product of dimensions on the
union of the operands' incident edges, counting the contracted edge once.
Equivalently, it is the two input element counts multiplied together and divided
by the shared bond dimension. Counts use checked `std::uint64_t` arithmetic;
fused dimensions must also fit the shape's dimension type.

These are scalar MACs, not real/complex FLOPs or measured time. Fusion changes
axis representation without adding dense MACs. Transpose traffic, copies,
communication, backend workspaces, and the lifetime of multiple tensors are not
modeled. The maximum of `step_elements` is the largest result tensor, not peak
allocated memory. The input graph must be simple, as required by the numerical
API; no outer-product operation is introduced.
