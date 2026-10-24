"""Application ownership policy: initial state, portable direct reads and logs."""
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

mpi, estimator, restart, complex_sampler, real_sampler, evidence_root = sys.argv[1:]
Path(evidence_root).mkdir(parents=True, exist_ok=True)
root = Path(tempfile.mkdtemp(prefix="accepted-", dir=evidence_root))
env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
count = 0

def run(binary, *args, ranks=2, ok=True, contains=None):
    global count
    count += 1
    command = [mpi, "--oversubscribe", "-np", str(ranks), binary, *map(str,args)]
    result = subprocess.run(command, env=env, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=60)
    (root / f"case-{count:02d}.log").write_text("COMMAND " + " ".join(command) + "\n" + result.stdout)
    assert (result.returncode == 0) == ok, result.stdout
    if contains: assert contains in result.stdout, result.stdout
    print(f"CASE {count} {'PASS' if ok else 'EXPECTED_REJECTION'} {' '.join(map(str,args))}")
    return result.stdout

def layout(output, source, expected):
    assert output.count("# site ownership: ") == 1, output
    assert f"# site ownership: {source}\n" in output
    rows = re.findall(r"^# MPI rank (\d+):([^\n]*)$", output, re.M)
    assert [[int(x) for x in sites.split()] for rank,sites in rows] == expected, rows

circuit = root / "prepare.qasm"
circuit.write_text('OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[4];\nx q[0];\ncx q[0],q[1];\ncx q[1],q[2];\ncx q[2],q[3];\nbarrier q;\n')
pauli = root / "pauli.txt"; pauli.write_text("1.0 0.0 IIIZ\n")
edges = root / "chain.edges"; edges.write_text("0 1\n1 2\n2 3\n")
lines = root / "chain.lines"; lines.write_text("0\n1\n2\n3\n")
owners = root / "owners.txt"; owners.write_text("# noncontiguous, deliberately uneven and unsorted\n3 0\n1 1 # comment\n0 0\n\n2 0\n")
other = root / "other.txt"; other.write_text("0 1\n1 0\n2 1\n3 1\n")
common = ["--circuit",circuit,"--sparse_pauli",pauli]
sampling = ["--linefile",lines,"--num_samples","16","--sampling_method","1","--bmps_bond_dim","16"]
for binary in (estimator,restart):
    extra = ["--measurement_barrier","0"] if binary == restart else []
    for path,expected in ((None,[[0,1],[2,3]]),(owners,[[0,2,3],[1]])):
        output = run(binary,*common,*extra,*(["--ownerfile",path] if path else []))
        layout(output,str(path) if path else "automatic",expected)
        values = re.findall(r"^ Expectation value of single site operator.* ([^ ]+)$",output,re.M)
        assert len(values)==1 and abs(float(values[0])+1)<1e-10,output

# Use a full-support complex state for sampling; zero-weight initial MCMC
# configurations need not escape under single-site updates.
circuit.write_text('OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[4];\nh q[0];\nh q[1];\nh q[2];\nh q[3];\ncz q[0],q[1];\ncz q[1],q[2];\ncz q[2],q[3];\ns q[3];\nbarrier q;\n')
saved = root / "saved"
run(restart,*common,"--ownerfile",owners,"--step_end","1","--portable_save",saved,
    "--savename",root / "legacy-")
# Actual tensors are read by the newly selected owners; roundtrip keeps every V/E byte.
copy = root / "copy"
output = run(restart,*common,"--ownerfile",other,"--portable_load",saved,"--step_end","1","--portable_save",copy)
layout(output,str(other),[[1],[0,2,3]])
assert "PORTABLE_READ rank=0 V=1" in output and "PORTABLE_READ rank=1 V=3" in output
for p in list(saved.glob("vertices/*"))+list(saved.glob("messages/*")):
    assert p.read_bytes()==(copy/p.relative_to(saved)).read_bytes()
output = run(complex_sampler,"--portable_input",saved,"--ownerfile",other,*sampling,"--num_samples","256")
layout(output,str(other),[[1],[0,2,3]])
samples = re.findall(r"^ sample \d+\.0: (\d+)$",output,re.M)
assert len(samples)==256
chi2=sum((samples.count(f"{i:04b}")-16)**2/16 for i in range(16))
assert chi2<70,(chi2,samples)
output = run(complex_sampler,"--portable_input",saved,*sampling)
layout(output,"automatic",[[0,1],[2,3]])
# The checkpoint was written with TWO ranks. Read it with FOUR different owners.
four=root/"four-ranks.txt";four.write_text("0 3\n1 1\n2 0\n3 2\n")
fourcopy=root/"four-rank-copy"
output=run(restart,*common,"--portable_load",saved,"--ownerfile",four,
           "--step_end","1","--portable_save",fourcopy,ranks=4)
layout(output,str(four),[[2],[1],[3],[0]])
for p in list(saved.glob("vertices/*"))+list(saved.glob("messages/*")):
    assert p.read_bytes()==(fourcopy/p.relative_to(saved)).read_bytes()
output=run(complex_sampler,"--portable_input",saved,"--ownerfile",four,
           *sampling,"--num_samples","256",ranks=4)
layout(output,str(four),[[2],[1],[3],[0]])
samples=re.findall(r"^ sample \d+\.0: (\d+)$",output,re.M)
assert len(samples)==256
assert sum((samples.count(f"{i:04b}")-16)**2/16 for i in range(16))<70
# Both random-state scalar types initialize directly on the requested owners.
for binary in (real_sampler,complex_sampler):
    output=run(binary,"--inputedge",edges,"--ownerfile",owners,"--p_bond_dim","1",*sampling)
    layout(output,str(owners),[[0,2,3],[1]])
    assert re.findall(r"^ sample \d+\.0: (\d+)$",output,re.M)==["0000"]*16

# Legacy inputs already carry a partition: matching maps are accepted; remaps are rejected.
for binary,flags in ((restart,[*common,"--loadname",root / "legacy-","--step_start","1","--step_end","1"]),
                     (complex_sampler,["--inputname",root / "legacy-","--inputedge",edges,*sampling])):
    output=run(binary,*flags,"--ownerfile",owners)
    layout(output,str(owners),[[0,2,3],[1]])
    output=run(binary,*flags)
    layout(output,"checkpoint",[[0,2,3],[1]])
    run(binary,*flags,"--ownerfile",other,ok=False,contains="differs from legacy checkpoint ownership")

# All apps reject malformed or incomplete files collectively before tensor work.
for kind,text,message in (
    ("duplicate","0 0\n0 1\n1 0\n2 0\n3 1\n","duplicate site"),
    ("missing","0 0\n1 1\n2 0\n","missing site"),
    ("unknown","0 0\n1 1\n2 0\n3 1\n4 0\n","unknown site"),
    ("range","0 0\n1 2\n2 0\n3 1\n","rank out of range"),
    ("negative","0 0\n1 -1\n2 0\n3 1\n","rank out of range"),
    ("extra","0 0 1\n","requires site_id mpi_rank"),
    ("number","0 nope\n","requires integer"),
):
    bad=root/f"{kind}.txt";bad.write_text(text)
    for binary,flags in ((estimator,common),(restart,common),(real_sampler,["--inputedge",edges,*sampling])):
        run(binary,*flags,"--ownerfile",bad,ok=False,contains=message)
for binary,flags in ((estimator,common),(restart,common),(real_sampler,["--inputedge",edges,*sampling])):
    run(binary,*flags,"--ownerfile",root/"absent.txt",ok=False,contains="cannot open ownerfile")
    run(binary,"--ownerfile",ok=False,contains="requires a value")
empty=root/"empty-rank.txt";empty.write_text("0 0\n1 0\n2 0\n3 0\n")
for binary,flags in ((estimator,common),(restart,common),
                     (real_sampler,["--inputedge",edges,*sampling])):
    output=run(binary,*flags,"--ownerfile",empty)
    layout(output,str(empty),[[0,1,2,3],[]])
run(complex_sampler,"--portable_input",saved,"--ownerfile",empty,*sampling,
    ok=False,contains="does not support empty ranks")
# Default partition uses the same get_range convention, including remainders.
six=root/"six.edges";six.write_text("0 1\n1 2\n2 3\n3 4\n4 5\n")
sixlines=root/"six.lines";sixlines.write_text("0\n1\n2\n3\n4\n5\n")
output=run(real_sampler,"--inputedge",six,"--linefile",sixlines,"--num_samples","1",ranks=4)
layout(output,"automatic",[[0,1],[2,3],[4],[5]])
print(f"PASS {count} ownership application cases; EVIDENCE {root}")
