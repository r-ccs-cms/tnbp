# Symbolic graph contraction costs

`graph_contraction_cost_symbolic` and `graph_schedule_cost_symbolic` accept an
ordered graph and `map<Site, vector<string>>` shapes. Each vector follows the
site's incident-edge order. `D` denotes an original virtual dimension and `x`
denotes chi: `DDx` means D² chi, and `1` is a dimension-one axis. Empty shapes
represent scalars. Dimension strings are products, not arbitrary expressions.

The result is a tuple of total MAC polynomial, per-step MAC powers, per-step
output-size powers, final edges, final symbolic shapes and per-step shared graph
operations (including common-neighbor fusion). All inputs are const;
invalid input and overflow leave them unchanged. The explicit schedule contains
`{site_a, site_b, site_c}` triples. It is not the edge-priority schedule accepted
by `graph_tensor_network_contraction_cost`.

The numerical and symbolic evaluators share the graph transformation, including
common-neighbor fusion and edge order. String parsing and initial shape checks
occur once per public evaluation. `graph_symbolic_cost_value` evaluates the
polynomial at positive integer model dimensions with checked arithmetic.
These costs count dense MACs, not data movement or peak allocated memory.

`build_graph_symbolic_region_plans` enumerates connected binary contractions of
up to five leaves in the full graph, retaining the target label and owner. It
returns prepared candidate-weighted cost polynomials plus a separate validated
baseline. The footprint includes external common-neighbor owners. Larger
regions use the baseline. `select_graph_symbolic_region_plan` compares costs at
known model dimensions, then maximum output size and lexicographic steps.
Unknown model dimensions or arithmetic overflow select the baseline. These are
model costs; they do not promise optimality for nonuniform truncated dimensions.

The default full-state heat-bath updater uses these plans on quotient graphs
of the original three-line edges. Original labels, edge order and native
common-neighbor fusion determine each step; no artificial external leaves are
introduced. The initial model uses the maximum input virtual dimension for D
and the configured boundary bond limit for chi. `mcmc_boundary_edge_model` consumes the ordered graph from `mcmc_init_edges`,
which already replays boundary absorption and bond shifts. Internal bonds of
absorbed/truncated boundaries use chi. The outermost unabsorbed row and edges
touching a target retain D. Actual truncated dimensions may differ. Standalone
line planners without checkpoint metadata assume supplied boundaries use chi;
they can receive explicit per-edge dimension strings instead.

At sampling entry, original input physical dimensions and model D are shared
once. All ranks deterministically prepare both sweep directions, environment
plans, window plans and the owner-set registry. The plans and communicators are
reused across all samples/sweeps and freed at return. Numerical region calls use
only their participating owners, including exterior common neighbors. There is
no line-wide shape reduction. A direct line call without prepared input metadata
uses the validated baseline and shares candidate counts within each window.

Within a region, `graph_tensor_contraction` performs the entire native step,
including exterior axis changes. Common subtrees run once; target-dependent
steps run per state. Private intermediate tensors pass directly to the next
step. Only retained environment/candidate cuts are normalized to unfused
original-edge order. This lets subsequent regions use a different valid native
fusion layout without recomputing environments. Factor dimensions needed to
normalize these cuts travel only between operand owners, once per step, and
are reused across candidates. Raw slices and boundary inputs remain immutable.

One-site windows retain A_s; two-site windows retain A_s and B_t and evaluate
all p*q scalar joins. Join metadata is compiled once per window. The selected
A_s is moved into the next prefix; the proper suffix is built once and consumed
as the sweep advances. No candidate reconstructs the exterior network. A
branch owns its mutable operands and common neighbors; unaffected exterior
tensors are never copied. Shared values are copied before destructive use and
moved on their final branch use. The previous prefix is consumed by A preparation;
the current suffix is consumed after its exterior-axis use, and B_t is consumed
at its last scalar join. These decisions include subsequent windows, while A_s
remains protected until roulette finishes. The graph API's own B-operand copy remains.

Numerical failures print a diagnostic and abort the MPI job. Scalar/weight
failures are notified once at window end; the selected state is broadcast only
within that window. Legacy sampling methods are unchanged. The implementation
uses free functions, tuples and containers, with no new owning class.

Region search parses and validates the shared string shapes once, then evaluates
all candidate schedules and the selected common-first order from that validated
input. Each evaluation owns its mutable graph/shape state; caller inputs remain
unchanged even on failure. Incident-axis validation scans edges in their stored
order without rescanning the full graph for each endpoint.

Scalar joins prepare owner mapping, local operand IDs and packed cut shapes once
per window, reusing them for every state pair. Canonically ordered single-edge
cuts require only reshape, with no per-pair factor map or cut normalization.
The mutable graph API still receives private working containers, and its normal
operand-owner communication and checks are unchanged.
