"""Small-system CLI acceptance, using only Python's standard library.
The native-record reader is intentionally CPU-specific; the application is not.
"""
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import struct
import subprocess
import sys
import tempfile

mpi, binary, evidence = sys.argv[1:4]
Path(evidence).mkdir(parents=True, exist_ok=True)
root = Path(tempfile.mkdtemp(prefix="accepted-", dir=evidence))
env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", TCAPI_VERBOSE="0")
checks = 0
runs = 0
metrics = {}

def check(ok, message):
    global checks
    if not ok:
        raise AssertionError(message)
    checks += 1

def put(name, text):
    path = root / name
    path.write_text(text)
    return path

def run(*args, ranks=1, ok=True):
    global runs
    runs += 1
    command = [mpi, "--oversubscribe", "-n", str(ranks), binary, *map(str, args)]
    result = subprocess.run(command, env=env, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=60)
    (root / f"run-{runs:03}.log").write_text(shlex.join(command) + "\n" + result.stdout)
    check((result.returncode == 0) == ok, f"run {runs}: unexpected exit\n{result.stdout}")
    if ok:
        for line in result.stdout.splitlines():
            check(bool(re.match(r"^(# )?\d{18} ", line)), f"missing timestamp: {line}")
        for line in result.stdout.splitlines():
            if re.match(r"^# ", line):
                check(bool(re.match(r"^# \d{18} ", line)), "input comment format")
    return result.stdout

def summaries(text):
    return [dict(re.findall(r"(\w+)=([^ ]+)", line)) for line in text.splitlines() if " SUMMARY " in line]

def tensor(path):
    raw = path.read_bytes()
    rank, = struct.unpack_from("=Q", raw)
    shape = struct.unpack_from("=" + "i"*rank, raw, 8)
    pos = 8+4*rank
    re_, im = struct.unpack_from("=dd", raw, pos)
    scale = complex(re_, im)
    pos += 16
    count = math.prod(shape)
    data = [complex(*struct.unpack_from("=dd", raw, pos+16*i))*scale for i in range(count)]
    check(pos+16*count == len(raw), "native record length")
    return shape, data

def normalize(v):
    norm = math.sqrt(sum(abs(x)**2 for x in v))
    return [x/norm for x in v]

def state(directory, n, edges):
    tensors = [tensor(directory / "vertices" / f"site-{i:06}.ten") for i in range(n)]
    dims = []
    for a, b in edges:
        incident = [e for e, ends in enumerate(edges) if a in ends]
        dims.append(tensors[a][0][incident.index(len(dims))])
    result = [0j]*(1 << n)
    for coordinates in itertools.product(*(range(d) for d in dims)):
        for bits in range(1 << n):
            value = 1+0j
            for site in range(n):
                shape, data = tensors[site]
                indices = [coordinates[e] for e, ends in enumerate(edges) if site in ends] + [(bits >> site) & 1]
                index, stride = 0, 1
                for c, d in zip(indices, shape):
                    index += c*stride
                    stride *= d
                value *= data[index]
            result[bits] += value
    return normalize(result)

def state_distance(a, b):
    overlap = sum(x.conjugate()*y for x, y in zip(a, b))
    phase = overlap / abs(overlap) if abs(overlap) else 1
    return math.sqrt(sum(abs(y-phase*x)**2 for x, y in zip(a, b)))

def same_state(a, b, tolerance=2e-10):
    return state_distance(a, b) < tolerance

def pauli(v, label):
    result = [0j]*len(v)
    for bits, amplitude in enumerate(v):
        out = bits
        for q, op in enumerate(reversed(label)):
            bit = (bits >> q) & 1
            if op in "XY":
                out ^= 1 << q
            if op == "Y":
                amplitude *= -1j if bit else 1j
            if op == "Z" and bit:
                amplitude = -amplitude
        result[out] += amplitude
    return result

def evolve(v, layers, schedule, order=2):
    actions = [(i, 1) for i in range(len(layers))]
    if order == 2 and layers:
        actions = [(i, .5) for i in range(len(layers)-1)] + [(len(layers)-1, 1)] + [(i, .5) for i in reversed(range(len(layers)-1))]
    for dt, count in schedule:
        for _ in range(count):
            for layer, fraction in actions:
                for coefficient, label in layers[layer]:
                    p = pauli(v, label)
                    x = dt*fraction*coefficient
                    v = [math.cosh(x)*a-math.sinh(x)*b for a, b in zip(v, p)]
                v = normalize(v)
    return v

def energy(v, terms):
    return sum(c*sum(a.conjugate()*b for a, b in zip(v, pauli(v, label))) for c, label in terms)

edges = put("edges.txt", "# reversed endpoints must be normalized\n1 0\n")
h = put("H.txt", ".23 0 IY\n-.37 0 XI\n.7 0 ZZ\n.11 0 II\n")
schedule = put("schedule.txt", ".08 3\n.025 4 # second stage\n")
base = ["--inputedge", edges, "--sparse_pauli", h, "--schedule", schedule, "--initial_state", "plus",
        "--max_bond_dim", "8", "--sv_min", "1e-12", "--truncation_error", "1e-12", "--bp_tolerance", "1e-12"]
terms = [(.23, "IY"), (-.37, "XI"), (.7, "ZZ"), (.11, "II")]
layers = [terms[:2], [terms[2]]]
full = root / "full"
text = run(*base, "--portable_save", full)
check("EDGES 1\n0 0 1\n" in (full/"manifest.txt").read_text(), "saved endpoints normalized")
rows = summaries(text)
check([int(x["step"]) for x in rows] == list(range(1, 8)), "default measurement every step")
check(all(x["bond_dim_max"] != "N/A" and x["trunc_error_max"] != "N/A" for x in rows), "graph statistics")
reference = evolve([.5]*4, layers, [(.08, 3), (.025, 4)])
metrics["two_site_exact_state_distance"] = state_distance(state(full, 2, [(0, 1)]), reference)
check(metrics["two_site_exact_state_distance"] < 2e-10, "full CLI state versus independent exact Trotter")
check(abs(float(rows[-1]["energy"])-energy(reference, terms).real) < 2e-10, "coefficient-weighted energy including Y and identity")
check(abs(float(rows[-1]["tau"])-.34) < 1e-15, "cumulative time")

split = root / "split"
run(*base, "--step_end", "3", "--measure_every", "2", "--portable_save", split)
resumed = root / "resumed"
text = run("--portable_load", split, "--step_end", "7", "--measure_every", "3", "--portable_save", resumed)
check([int(x["step"]) for x in summaries(text)] == [6, 7], "absolute skip cadence and terminal measurement")
check(same_state(state(resumed, 2, [(0, 1)]), reference), "same-rank restart state")
check((full / "COMMITTED").read_bytes() == (resumed / "COMMITTED").read_bytes(), "same-rank restart state digest identical")
for a in (full / "vertices").glob("*.ten"):
    check(a.read_bytes() == (resumed / "vertices" / a.name).read_bytes(), "same-rank V payload identical")
for a in (full / "messages").glob("*.ten"):
    check(a.read_bytes() == (resumed / "messages" / a.name).read_bytes(), "same-rank E payload identical")

owners = put("owners.txt", "0 1\n1 0\n")
parallel = root / "parallel"
text = run("--portable_load", split, "--ownerfile", owners, "--portable_save", parallel, ranks=2)
check(same_state(state(parallel, 2, [(0, 1)]), reference), "changed-rank custom-owner continuation")
check(abs(float(summaries(text)[-1]["energy"])-energy(reference, terms).real) < 2e-10, "changed-rank energy")

skip = root / "skip"
text = run(*base, "--measure_every", "1000", "--portable_save", skip)
check([int(x["step"]) for x in summaries(text)] == [7], "measurement skip with terminal output")
check((full / "COMMITTED").read_bytes() == (skip / "COMMITTED").read_bytes(), "measurement cadence does not change V/E")
check(summaries(text)[0]["window_start"] == "0" and summaries(text)[0]["window_end"] == "7", "statistics window")
skip_row = summaries(text)[0]
check(float(skip_row["trunc_error_min"]) == min(float(r["trunc_error_min"]) for r in rows), "interval error minimum")
check(float(skip_row["trunc_error_max"]) == max(float(r["trunc_error_max"]) for r in rows), "interval error maximum")
check(skip_row["bond_dim_min"] == rows[-1]["bond_dim_min"] and skip_row["bond_dim_max"] == rows[-1]["bond_dim_max"], "bond dimensions describe current state")

# Restart inside a schedule stage as well as at its boundary.
inside = root / "inside"
run(*base, "--step_end", "2", "--portable_save", inside)
inside_final = root / "inside-final"
run("--portable_load", inside, "--portable_save", inside_final)
check((full / "COMMITTED").read_bytes() == (inside_final / "COMMITTED").read_bytes(), "restart within a stage retains schedule position")

# An explicit first-order calculation and a two-site Y measurement.
h_two_y = put("H-two-y.txt", ".4 0 YZ\n-.2 0 IX\n")
first_order = root / "first-order"
text = run("--inputedge", edges, "--sparse_pauli", h_two_y, "--time_step", ".02", "--num_steps", "4",
           "--trotter_order", "1", "--sv_min", "1e-12", "--truncation_error", "1e-12", "--portable_save", first_order)
expected_first = evolve([1,0,0,0], [[(-.2,"IX")],[(.4,"YZ")]], [(.02,4)], order=1)
check(same_state(state(first_order,2,[(0,1)]),expected_first), "first-order evolution from default zero state")
check(abs(float(summaries(text)[-1]["energy"])-energy(expected_first,[(.4,"YZ"),(-.2,"IX")]).real)<2e-10, "two-site Y measurement")


no_op = root / "no-op"
text = run("--portable_load", split, "--step_end", "3", "--portable_save", no_op)
check(len(summaries(text)) == 1 and summaries(text)[0]["trunc_error_max"] == "N/A", "no evolution summary")
check((split / "COMMITTED").read_bytes() == (no_op / "COMMITTED").read_bytes(), "no evolution payload unchanged")

# State-only input and explicit warm start reset the application's clock.
state_only = root / "state-only"
shutil.copytree(split, state_only)
(state_only / "imaginary_time.resume").unlink()
new = ["--sparse_pauli", h, "--time_step", ".01", "--num_steps", "2"]
warm = root / "warm"
text = run("--portable_load", split, "--restart_mode", "initial", *new, "--portable_save", warm)
check([int(x["step"]) for x in summaries(text)] == [1, 2], "warm start clock reset")
plain = root / "plain"
run("--portable_load", state_only, *new, "--portable_save", plain)
check((warm / "COMMITTED").read_bytes() == (plain / "COMMITTED").read_bytes(), "state-only and warm start agree")

# Isolated sites, basis-file initialization, one-site evolution and empty edges.
empty = put("empty.txt", "# no edges\n")
h1 = put("H1.txt", ".6 0 Y\n-.2 0 Z\n.13 0 I\n")
bits = put("bits.txt", "# site bit\n0 1\n")
iso = root / "iso"
text = run("--inputedge", empty, "--sparse_pauli", h1, "--initial_statefile", bits,
           "--time_step", ".03", "--num_steps", "5", "--portable_save", iso)
expected = evolve([0, 1], [[(.6, "Y")], [(-.2, "Z")]], [(.03, 5)])
check(same_state(state(iso, 1, []), expected), "isolated one-site evolution")
check(summaries(text)[-1]["bond_dim_max"] == "N/A" and summaries(text)[-1]["trunc_error_max"] == "N/A", "empty-edge statistics")
check(abs(float(summaries(text)[-1]["energy"])-energy(expected, [(.6,"Y"),(-.2,"Z"),(.13,"I")]).real) < 2e-12, "isolated measurement")
iso_resume = root / "iso-resume"
run("--portable_load", iso, "--portable_save", iso_resume)
check((iso / "COMMITTED").read_bytes() == (iso_resume / "COMMITTED").read_bytes(), "isolated portable roundtrip")

# One connected component plus an isolated qubit; a long-range term is evolved,
# but its energy must not be silently reported as the complete Hamiltonian.
e3 = put("e3.txt", "1 0\n1 2\n")
h4 = put("H4.txt", ".2 0 ZIIZ\n-.3 0 YIII\n.15 0 IIXI\n")
mixed = root / "mixed"
text = run("--inputedge", e3, "--sparse_pauli", h4, "--initial_state", "plus", "--time_step", ".04", "--num_steps", "3",
           "--sv_min", "1e-12", "--truncation_error", "1e-12", "--portable_save", mixed, ok=False)
# ZIIZ spans an isolated site and is intentionally rejected; use connected long range instead.
h4 = put("H4valid.txt", ".2 0 IZIZ\n-.3 0 YIII\n.15 0 IIXI\n")
text = run("--inputedge", e3, "--sparse_pauli", h4, "--initial_state", "plus", "--time_step", ".04", "--num_steps", "3",
           "--sv_min", "1e-12", "--truncation_error", "1e-12", "--portable_save", mixed)
check(summaries(text)[-1]["energy"] == "N/A" and summaries(text)[-1]["unsupported_terms"] == "1", "partial-energy label")
check("partial_energy" in summaries(text)[-1], "partial-energy field")
# All these terms commute. Leaf-side compression with BP between batches
# should match exact evolution in this tree example. The previous all-edge
# frozen-environment approximation remains documented separately; no claim of
# exactness for arbitrary loop graphs or truncation parameters is implied.
expected = evolve([.25]*16, [[(.15,"IIXI"),(.2,"IZIZ"),(-.3,"YIII")]], [(.04,3)])
mixed_state = state(mixed,4,[(0,1),(1,2)])
metrics["long_range_repeated_approximation_distance"] = state_distance(mixed_state,expected)
check(metrics["long_range_repeated_approximation_distance"]<2e-10, "leaf-batch BP long-range evolution versus exact state")
check(abs(energy(mixed_state,[(1,"YIII")])-energy(expected,[(1,"YIII")]))<2e-10, "independent isolated component remains exact")

mixed_first = root / "mixed-first"
run("--inputedge", e3, "--sparse_pauli", h4, "--initial_state", "plus", "--time_step", ".04", "--num_steps", "3",
    "--step_end", "1", "--sv_min", "1e-12", "--truncation_error", "1e-12", "--portable_save", mixed_first)
# Endpoint normalization is allowed at input, but edge IDs in a saved state
# must never be reassigned by sorting or deduplicating its graph.
equivalent_e3=put("equivalent-e3.txt","0 1\n1 2\n")
run("--portable_load",mixed_first,"--inputedge",equivalent_e3,"--step_end","1")
reordered_e3=put("reordered-e3.txt","1 2\n0 1\n")
text=run("--portable_load",mixed_first,"--inputedge",reordered_e3,"--step_end","1",ok=False)
check("graph/order differs from portable state" in text,"same edges in different order rejected on resume")
duplicate_e3=put("duplicate-e3.txt","0 1\n1 2\n1 0\n")
text=run("--inputedge",duplicate_e3,"--sparse_pauli",h4,"--time_step",".04","--num_steps","1",ok=False)
check("duplicate undirected edge" in text,"duplicate input rejected rather than silently uniqued")
expected_first = evolve([.25]*16, [[(.15,"IIXI"),(.2,"IZIZ"),(-.3,"YIII")]], [(.04,1)])
check(same_state(state(mixed_first,4,[(0,1),(1,2)]),expected_first,2e-10), "first long-range layer including isolated component versus exact operator")
mixed_restart = root / "mixed-restart"
run("--portable_load", mixed_first, "--portable_save", mixed_restart)
check((mixed / "COMMITTED").read_bytes() == (mixed_restart / "COMMITTED").read_bytes(), "long-range approximate evolution resumes identically")
mixed_parallel = root / "mixed-parallel"
run("--portable_load", mixed_first, "--portable_save", mixed_parallel, ranks=2)
check(same_state(mixed_state,state(mixed_parallel,4,[(0,1),(1,2)]),2e-9), "long-range MPI continuation preserves the existing approximation")

# Leaf-side compression on a branching tree, and on two independent trees.
for name,n,graph,terms_,expected_batches in [
    ("branch",4,[(0,1),(1,2),(1,3)],[(.17,"ZYIZ")],[[1],[2],[0]]),
    ("forest",7,[(0,1),(1,2),(1,3),(4,5),(5,6)],[(.17,"IIIZYIZ"),(-.11,"YZYIIII")],[[1,4],[2,3],[0]])]:
    graph_file=put(name+"-edges.txt","".join(f"{a} {b}\n" for a,b in graph))
    ham=put(name+"-H.txt","".join(f"{c} 0 {label}\n" for c,label in terms_))
    args=["--inputedge",graph_file,"--sparse_pauli",ham,"--initial_state","plus",
          "--time_step",".04","--num_steps","10","--sv_min","1e-12","--truncation_error","1e-12","--bp_tolerance","1e-12"]
    final=root/(name+"-final")
    text=run(*args,"--measure_every","3","--portable_save",final)
    found=[list(map(int,m.group(1).split())) for m in re.finditer(r"compression_batch layer=0 batch=\d+ edges=([^\n]+)",text)]
    check(found==expected_batches,name+" deterministic leaf-side batches")
    check("BP_NOT_CONVERGED" not in text,name+" all BP stages converge")
    exact=evolve([1/math.sqrt(1<<n)]*(1<<n),[terms_],[(.04,10)])
    psi=state(final,n,graph)
    metrics[name+"_leaf_batch_distance"]=state_distance(psi,exact)
    check(same_state(psi,exact),name+" leaf-batch/BP versus independent exact evolution")
    check(summaries(text)[-1]["energy"]=="N/A",name+" general Pauli measurement remains unsupported")
    prefix=root/(name+"-prefix")
    run(*args,"--step_end","4","--portable_save",prefix)
    resumed=root/(name+"-resumed")
    run("--portable_load",prefix,"--portable_save",resumed)
    check((final/"COMMITTED").read_bytes()==(resumed/"COMMITTED").read_bytes(),name+" same-rank restart preserves V/E")
    noop=root/(name+"-noop-two-ranks")
    run("--portable_load",prefix,"--step_end","4","--portable_save",noop,ranks=2)
    check((prefix/"COMMITTED").read_bytes()==(noop/"COMMITTED").read_bytes(),name+" changed-rank zero-step save preserves every record")
    parallel=root/(name+"-parallel")
    run("--portable_load",prefix,"--portable_save",parallel,ranks=2)
    check(same_state(state(parallel,n,graph),exact),name+" changed-rank continuation versus exact state")
    whole_scale=.04*10*sum(abs(c) for c,_ in terms_)
    for cp in [final,resumed,parallel]:
        value=float(re.search(r"OPERATOR_LOG_SCALE ([^\n]+)",(cp/"imaginary_time.resume").read_text()).group(1))
        check(abs(value-whole_scale)<2e-15,name+" complete operator scale independent of rank count and restart")

    if name=="branch":
        # A zero iteration budget exposes the position of every BP invocation
        # without adding production tracing or any extra convergence calls.
        text=run(*args,"--max_bp_iterations","0","--step_end","1")
        warned=[int(v) for v in re.findall(r"BP_NOT_CONVERGED step=1 layer=0 before_batch=(\d+)",text)]
        check(warned==[0,1,2],"BP occurs before first batch and between batches, never after last or for measurement")
        two_full=root/"branch-two-full"
        run(*args,"--portable_save",two_full,ranks=2)
        two_prefix=root/"branch-two-prefix"
        run(*args,"--step_end","4","--portable_save",two_prefix,ranks=2)
        two_resumed=root/"branch-two-resumed"
        run("--portable_load",two_prefix,"--portable_save",two_resumed,ranks=2)
        check((two_full/"COMMITTED").read_bytes()==(two_resumed/"COMMITTED").read_bytes(),"MPI2 continuous and resumed V/E identical")
        reverse=root/"branch-two-to-one"
        run("--portable_load",two_prefix,"--portable_save",reverse)
        check(same_state(state(reverse,n,graph),exact),"MPI2 to MPI1 continuation versus exact state")

# Two edge-disjoint strings share only the central identity routing site.
star_graph=[(0,1),(0,2),(0,3),(0,4)]
star_edges=put("star-edges.txt","".join(f"{a} {b}\n" for a,b in star_graph))
star_terms=[(.21,"IIYXI"),(-.16,"ZZIII")]
star_h=put("star-H.txt","".join(f"{c} 0 {label}\n" for c,label in star_terms))
star_args=["--inputedge",star_edges,"--sparse_pauli",star_h,"--initial_state","plus","--time_step",".03","--num_steps","8","--sv_min","1e-12","--truncation_error","1e-12","--bp_tolerance","1e-12"]
star_final=root/"star-final"
text=run(*star_args,"--portable_save",star_final)
check(" layers: 2" in text,"shared identity routing site splits CLI layers")
star_batches=[(int(m.group(1)),list(map(int,m.group(2).split()))) for m in re.finditer(r"compression_batch layer=(\d+) batch=\d+ edges=([^\n]+)",text)]
check(star_batches==[(0,[1]),(0,[0]),(1,[3]),(1,[2])],"shared-site strings stay separate and leaf-first")
star_exact=evolve([1/math.sqrt(32)]*32,[[star_terms[0]],[star_terms[1]]],[(.03,8)])
metrics["shared_routing_site_split_distance"]=state_distance(state(star_final,5,star_graph),star_exact)
check(metrics["shared_routing_site_split_distance"]<2e-10,"site-exclusive star evolution versus independent exact Trotter")
star_prefix=root/"star-prefix"
run(*star_args,"--step_end","3","--portable_save",star_prefix)
star_resumed=root/"star-resumed"
run("--portable_load",star_prefix,"--portable_save",star_resumed)
check((star_final/"COMMITTED").read_bytes()==(star_resumed/"COMMITTED").read_bytes(),"new site-disjoint plan restart preserves V/E")
star_parallel=root/"star-parallel"
run("--portable_load",star_prefix,"--portable_save",star_parallel,ranks=2)
check(same_state(state(star_parallel,5,star_graph),star_exact),"site-exclusive split continuation with changed MPI owners")

# Single-site layers on connected vertices still need local normalization.
onsite_h=put("onsite-connected.txt","-.3 0 IX\n.2 0 ZI\n")
onsite_cp=root/"onsite-connected"
text=run("--inputedge",edges,"--sparse_pauli",onsite_h,"--initial_state","plus","--time_step",".1","--num_steps","20","--portable_save",onsite_cp)
onsite_exact=evolve([.5]*4,[[(-.3,"IX"),(.2,"ZI")]],[(.1,20)])
check(same_state(state(onsite_cp,2,[(0,1)]),onsite_exact),"connected onsite-only evolution without edge compression")
check(summaries(text)[-1]["trunc_error_max"]=="N/A" and summaries(text)[-1]["bond_dim_max"]=="1","unused existing edges retain dimensions but have no truncation samples")

# A thousand cheap one-site steps test cadence and actual ground-state approach.
ground_h = put("ground.txt", "-1 0 Z\n")
ground_cp = root / "ground"
text = run("--inputedge",empty,"--sparse_pauli",ground_h,"--initial_state","plus","--time_step",".01","--num_steps","1000",
           "--measure_every","1000","--portable_save",ground_cp)
check(len(summaries(text))==1 and summaries(text)[0]["step"]=="1000", "1000-step measurement skip")
check(same_state(state(ground_cp,1,[]),[1,0],1e-8), "one-site ground-state approach")
check(abs(float(summaries(text)[0]["energy"])+1)<1e-12,"one-site ground energy")

identity_h = put("identity.txt", "1.2 0 II\n")
identity_cp = root / "identity"
text = run("--inputedge", edges, "--sparse_pauli", identity_h, "--time_step", ".1", "--num_steps", "3", "--portable_save", identity_cp)
check(same_state(state(identity_cp,2,[(0,1)]),[1,0,0,0]), "identity-only state unchanged")
check(float(summaries(text)[-1]["energy"]) == 1.2 and summaries(text)[-1]["trunc_error_max"] == "N/A", "identity-only energy and absent truncation")
check(abs(float(re.search(r"operator_log_scale=([^ ]+)", text).group(1))+.36) < 1e-14, "identity scalar once per step")

# Corruption, changed conditions and transactional publication errors.
run("--portable_load", split, "--time_step", ".02", "--num_steps", "2", ok=False)
run("--portable_load", split, "--max_bond_dim", "9", ok=False)
run("--portable_load", split, "--step_end", "2", ok=False)
run("--portable_load", split, "--step_end", "100", ok=False)
run("--portable_load", split, "--measure_every", "0", ok=False)
run(*base, "--portable_save", full, ok=False)
corrupt = root / "corrupt"
shutil.copytree(split,corrupt)
f=corrupt/"imaginary_time.resume"
f.write_text(f.read_text().replace("NEXT_STEP 3", "NEXT_STEP 2"))
run("--portable_load",corrupt,ok=False)
run("--portable_load",corrupt,"--restart_mode","initial",*new,ok=False)
pending = root / "pending"
shutil.copytree(split,pending)
(pending/"imaginary_time.resume.pending").write_text("incomplete")
run("--portable_load",pending,ok=False)
equivalent_graph=put("equivalent-edges.txt","0 1\n")
run("--portable_load",split,"--inputedge",equivalent_graph,"--step_end","3")
wrong_graph=put("wrong-edges.txt","")
run("--portable_load",split,"--inputedge",wrong_graph,ok=False)
bad_bits=put("bad-bits.txt","0 2\n1 0\n")
run("--inputedge",edges,"--sparse_pauli",h,"--time_step",".1","--num_steps","1","--initial_statefile",bad_bits,ok=False)
run("--inputedge",edges,"--sparse_pauli",h,"--time_step","0","--num_steps","1",ok=False)
# App metadata is independently integrity checked and bound to the state.
def modified_resume(name, old, new):
    target = root/name
    shutil.copytree(split,target)
    path = target/"imaginary_time.resume"
    text = path.read_text()
    body = text[:text.rfind("SHA256 ")]
    check(old in body,"metadata mutation target")
    body = body.replace(old,new)
    path.write_text(body+"SHA256 "+hashlib.sha256(body.encode()).hexdigest()+"\n")
    return target

run("--portable_load",modified_resume("old-layer-plan","PLAN_VERSION 2","PLAN_VERSION 1"),ok=False)
run("--portable_load",modified_resume("bad-compression","COMPRESSION leaf-batches-v1","COMPRESSION unknown"),ok=False)
run("--portable_load",modified_resume("old-draft","TNBP_IMAGINARY_TIME_RESUME 2","TNBP_IMAGINARY_TIME_RESUME 1"),ok=False)
run("--portable_load",modified_resume("bad-version","TNBP_IMAGINARY_TIME_RESUME 2","TNBP_IMAGINARY_TIME_RESUME 3"),ok=False)
run("--portable_load",modified_resume("bad-phase","PHASE step-complete","PHASE mid-step"),ok=False)
old_plan = re.search(r'PLAN "([0-9a-f]+)"',(split/"imaginary_time.resume").read_text()).group(1)
run("--portable_load",modified_resume("bad-plan",old_plan,"0"*64),ok=False)
wrong_state = root/"wrong-state"
shutil.copytree(full,wrong_state)
shutil.copy2(split/"imaginary_time.resume",wrong_state/"imaginary_time.resume")
run("--portable_load",wrong_state,ok=False)
staging_target=root/"staged"
Path(str(staging_target)+".imaginary-time.pending").mkdir()
run(*base,"--step_end","0","--portable_save",staging_target,ok=False)
check(not staging_target.exists(),"staging failure never publishes a final directory")
# Identical explicit conditions on continuation are accepted.
run("--portable_load",split,"--sparse_pauli",h,"--schedule",schedule,"--max_bond_dim","8","--step_end","3")
for key,value in [("--time_step","nan"),("--num_steps","-1"),("--trotter_order","3"),("--max_bond_dim","0")]:
    args=["--inputedge",edges,"--sparse_pauli",h,"--time_step",".1","--num_steps","1"]
    if key in args: args[args.index(key)+1]=value
    else: args += [key,value]
    run(*args,ok=False)
text=run("--inputedge",root,"--sparse_pauli",identity_h,"--time_step",".1","--num_steps","1",ok=False)
check("expected edge file, got directory" in text,"directory rejected as read error, not an unrelated Hamiltonian error")
for i,content in enumerate(["oops 1\n", "9999999999999999999999999 1\n"]):
    malformed=put(f"malformed-edge-{i}.txt",content)
    text=run("--inputedge",malformed,"--sparse_pauli",identity_h,"--time_step",".1","--num_steps","1",ok=False)
    check("invalid edge endpoint at line 1" in text,"malformed edge rejected before silently changing graph")
run(*base,"--time_step num_steps","1",ok=False)
run("--help")
(root/"metrics.json").write_text(json.dumps(metrics,indent=2)+"\n")
print("METRICS "+json.dumps(metrics,sort_keys=True))
print(f"PASS cli runs={runs} checks={checks} evidence={root}")
