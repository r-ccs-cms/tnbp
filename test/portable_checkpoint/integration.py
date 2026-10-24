"""CPU acceptance: real estimator output, changed MPI ownership, and corruption.
Uses only Python's standard library; no library is required by the product.
"""
import hashlib
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

mpi, fixture, estimator, sampler, real_sampler = sys.argv[1:6]
env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
with tempfile.TemporaryDirectory(prefix="tnbp-portable-") as tmp:
    root = Path(tmp)
    counter = 0
    transcripts = []

    def run(ranks, binary, *args, ok=True, contains=None):
        global counter
        counter += 1
        command = [mpi, "--oversubscribe", "-np", str(ranks), binary, *map(str, args)]
        result = subprocess.run(command, env=env, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=90)
        transcripts.append((counter, command, result.stdout))
        print(f"CASE {counter} ranks={ranks} {' '.join(map(str, args))} exit={result.returncode}")
        if (result.returncode == 0) != ok or (contains and contains not in result.stdout):
            print(result.stdout)
            raise AssertionError(command)
        return result.stdout

    def change_manifest(directory, transform):
        path = directory / "manifest.txt"
        path.write_text(transform(path.read_text()))
        (directory / "COMMITTED").write_text("TNBP_COMMITTED_V1 " + hashlib.sha256(path.read_bytes()).hexdigest() + "\n")

    def change_resume(directory, transform):
        path = directory / "estimator.resume"
        body = path.read_text().split("SHA256 ")[0]
        body = transform(body)
        path.write_text(body + "SHA256 " + hashlib.sha256(body.encode()).hexdigest() + "\n")

    for writers in (1, 4):
        original = root / f"fixture{writers}"
        run(writers, fixture, "save", original)
        for readers in (1, 2, 3, 4):
            run(readers, fixture, "load", original, root / f"roundtrip{writers}-{readers}", contains="VERIFIED_BYTES")
        run(4, fixture, "reverse", original, contains="VERIFIED_BYTES")
    run(2, fixture, "mismatch", root / "mismatch", ok=False, contains="boundary copies disagree")
    run(2, fixture, "limit", root / "fixture1", ok=False, contains="allocation limit")
    run(2, fixture, "bad-owner", root / "fixture1", ok=False, contains="owner coverage")
    run(5, fixture, "load", root / "fixture1", ok=False, contains="empty ranks")
    run(1, fixture, "save", root / "fixture1", ok=False, contains="already exists")
    for kind in ("corrupt", "missing", "pending", "shape", "duplicate", "codec", "manifest"):
        directory = root / kind
        shutil.copytree(root / "fixture1", directory)
        payload = directory / "vertices/site--00004.ten"
        if kind == "manifest": (directory / "manifest.txt").write_text("invalid")
        elif kind == "corrupt":
            data = bytearray(payload.read_bytes()); data[-1] ^= 1; payload.write_bytes(data)
        elif kind == "missing": payload.unlink()
        elif kind == "pending": (directory / "COMMITTED").unlink()
        elif kind == "shape":
            change_manifest(directory, lambda s: s.replace('"vertices/site--00004.ten" 104', '"vertices/site--00004.ten" 103'))
            # Always invalidate an indexed shape, independently of codec byte sizes.
            change_manifest(directory, lambda s: re.sub(r'("vertices/site--00004.ten" \d+ [a-f0-9]{64} 2) 2 2', r'\1 3 2', s))
        elif kind == "duplicate":
            change_manifest(directory, lambda s: s.replace('"vertices/site-000002.ten"', '"vertices/site--00004.ten"'))
        elif kind == "codec": change_manifest(directory, lambda s: s.replace("gqten-native-v1", "unsupported"))
        run(2, fixture, "load", directory, ok=False)

    # A four-site cluster state with complex phases and full sampling support.
    circuit = root / "cluster.qasm"
    circuit.write_text('OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[4];\nh q[0];\nh q[1];\nh q[2];\nh q[3];\ncz q[0],q[1];\ncz q[1],q[2];\ncz q[2],q[3];\ns q[3];\nbarrier q;\n')
    pauli = root / "pauli.txt"
    pauli.write_text("1.0 0.0 IIIZ\n")
    lines = root / "lines.txt"
    lines.write_text("0\n1\n2\n3\n")
    common = ["--circuit", circuit, "--sparse_pauli", pauli, "--max_bond_dim", "16"]
    sampler_reference = {}
    for writers in (1, 2, 4):
        saved = root / f"estimator{writers}"
        output = run(writers, estimator, *common, "--step_end", "1", "--portable_save", saved, contains="PORTABLE_SAVED")
        digest = re.search(r"PORTABLE_SAVED digest=(\w+)", output)[1]
        for readers in (1, 2, 4):
            output = run(readers, sampler, "--portable_input", saved, "--linefile", lines,
                         "--num_samples", "256", "--sampling_method", "1", "--bmps_bond_dim", "16", contains="PORTABLE_INPUT")
            assert f"digest={digest}" in output
            samples = re.findall(r"^ sample \d+\.0: (\d+)$", output, re.M)
            assert len(samples) == 256 and all(len(x) == 4 for x in samples)
            if writers == 1: sampler_reference[readers] = samples
            # Exact cluster-state computational-basis distribution is uniform.
            chi2 = sum((samples.count(f"{i:04b}") - 16) ** 2 / 16 for i in range(16))
            assert chi2 < 70, (chi2, samples)
            print(f"DISTRIBUTION writers={writers} readers={readers} n=256 chi2={chi2:.3f} digest={digest}")
        # Loading without evolution preserves all actual V/E bytes on a new partition.
        copied = root / f"estimator-roundtrip{writers}"
        output = run(4 if writers == 1 else 1, estimator, *common,
                     "--portable_load", saved, "--step_end", "1", "--portable_save", copied)
        assert f"digest={digest}" in output
        for payload in list(saved.glob("vertices/*")) + list(saved.glob("messages/*")):
            assert payload.read_bytes() == (copied / payload.relative_to(saved)).read_bytes()
    saved = root / "estimator1"
    # No E files exist: sampler must not stat/open/hash/decode them.
    vertices_only = root / "vertices-only"
    shutil.copytree(saved, vertices_only)
    shutil.rmtree(vertices_only / "messages")
    for readers in (1, 2, 4):
        run(readers, fixture, "vertices-only", vertices_only, contains="VERIFIED_BYTES")
        output = run(readers, sampler, "--portable_input", vertices_only, "--linefile", lines,
                     "--num_samples", "256", "--sampling_method", "1", "--bmps_bond_dim", "16")
        samples = re.findall(r"^ sample \d+\.0: (\d+)$", output, re.M)
        assert samples == sampler_reference[readers]
    run(2, estimator, *common, "--portable_load", vertices_only, "--step_end", "1", ok=False, contains="file length mismatch")
    corrupt_e = root / "corrupt-e"
    shutil.copytree(saved, corrupt_e)
    edge_file = next((corrupt_e / "messages").glob("*.ten"))
    data = bytearray(edge_file.read_bytes()); data[-1] ^= 1; edge_file.write_bytes(data)
    run(2, fixture, "vertices-only", corrupt_e, contains="VERIFIED_BYTES")
    run(2, sampler, "--portable_input", corrupt_e, "--linefile", lines, "--num_samples", "4", "--sampling_method", "1", "--bmps_bond_dim", "16")
    run(2, estimator, *common, "--portable_load", corrupt_e, "--step_end", "1", ok=False, contains="checksum mismatch")
    run(2, real_sampler, "--portable_input", saved, "--linefile", lines, ok=False, contains="backend/dtype/ABI mismatch")
    external = root / "wrong.edges"
    external.write_text("1 2\n0 1\n2 3\n")
    run(2, sampler, "--portable_input", saved, "--inputedge", external, "--linefile", lines, ok=False, contains="external graph differs")
    run(2, estimator, *common, "--portable_load", saved, "--step_end", "2", "--sv_min", "1e-7", ok=False, contains="resume conditions differ")
    unknown_app = root / "unknown-app"; shutil.copytree(saved, unknown_app)
    change_resume(unknown_app, lambda s: s.replace("TNBP_ESTIMATOR_RESUME 1", "TNBP_ESTIMATOR_RESUME 999"))
    run(2, sampler, "--portable_input", unknown_app, "--linefile", lines, "--num_samples", "4", "--sampling_method", "1", "--bmps_bond_dim", "16")
    run(2, estimator, *common, "--portable_load", unknown_app, "--step_end", "2", ok=False, contains="unsupported or malformed resume file")
    run(2, estimator, *common, "--portable_load", saved, "--step_start", "0", "--step_end", "2", ok=False, contains="step_start conflicts")
    no_app = root / "no-app"; shutil.copytree(saved, no_app)
    (no_app / "estimator.resume").unlink()
    run(2, sampler, "--portable_input", no_app, "--linefile", lines, "--num_samples", "4", "--sampling_method", "1", "--bmps_bond_dim", "16")
    run(2, estimator, *common, "--portable_load", no_app, "--step_end", "1",
        contains="PORTABLE_RESTART circuit=unspecified step_start=0")
    # Genuine next-cycle evolution; input is immutable and its manifest remains valid.
    run(2, estimator, *common, "--portable_load", saved, "--step_end", "2", "--portable_save", root / "evolved",
        contains="PORTABLE_RESTART circuit=same step_start=1")
    # Changed circuit starts a new step sequence on the saved graph and state.
    inverse = root / "inverse.qasm"
    inverse.write_text('OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[4];\nsdg q[3];\ncz q[2],q[3];\ncz q[1],q[2];\ncz q[0],q[1];\nh q[3];\nh q[2];\nh q[1];\nh q[0];\nbarrier q;\n')
    reset = root / "inverse-state"
    output = run(4, estimator, "--circuit", inverse, "--sparse_pauli", pauli,
                 "--portable_load", saved, "--step_start", "0", "--step_end", "1",
                 "--measurement_barrier", "0", "--portable_save", reset,
                 contains="PORTABLE_RESTART circuit=changed step_start=0")
    # U^-1 U |0000> must be zero: check actual state, not just successful loading.
    output = run(2, sampler, "--portable_input", reset, "--linefile", lines,
                 "--num_samples", "32", "--sampling_method", "1", "--bmps_bond_dim", "16")
    assert re.findall(r"^ sample \d+\.0: (\d+)$", output, re.M) == ["0000"] * 32
    assert "NEXT_STEP 1\n" in (reset / "estimator.resume").read_text()
    assert 'CIRCUIT "' + hashlib.sha256(inverse.read_bytes()).hexdigest() + '"' in (reset / "estimator.resume").read_text()
    run(2, estimator, "--circuit", inverse, "--sparse_pauli", pauli,
        "--portable_load", saved, "--step_start", "1", "--step_end", "2",
        ok=False, contains="step_start conflicts")
    run(2, estimator, "--circuit", inverse, "--sparse_pauli", pauli,
        "--portable_load", reset, "--step_end", "1", contains="PORTABLE_RESTART circuit=same step_start=1")
    # A circuit with no two-site gates must still use all saved edges and tensors.
    flip = root / "flip.qasm"
    flip.write_text('OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[4];\nx q[0];\nbarrier q;\n')
    flipped = root / "flipped"
    run(2, estimator, "--circuit", flip, "--sparse_pauli", pauli,
        "--portable_load", reset, "--step_end", "1", "--portable_save", flipped,
        contains="PORTABLE_RESTART circuit=changed step_start=0")
    output = run(4, sampler, "--portable_input", flipped, "--linefile", lines,
                 "--num_samples", "32", "--sampling_method", "1", "--bmps_bond_dim", "16")
    assert re.findall(r"^ sample \d+\.0: (\d+)$", output, re.M) == ["1000"] * 32
    bad = root / "outside.qasm"
    bad.write_text(flip.read_text().replace("q[4]", "q[5]").replace("x q[0]", "x q[4]"))
    run(2, estimator, "--circuit", bad, "--sparse_pauli", pauli,
        "--portable_load", saved, "--step_end", "1", ok=False,
        contains="circuit references site outside checkpoint graph")
    assert "NEXT_STEP 2\n" in (root / "evolved/estimator.resume").read_text()
    # Library-created state has no estimator metadata, and can initialize either app.
    for directory in (root / "fixture1", saved):
        manifest = (directory / "manifest.txt").read_text()
        assert manifest.startswith("TNBP_PORTABLE_TPS 3\n")
        assert "APPLICATION" not in manifest and "estimator" not in manifest and "WRITERS" not in manifest
    assert not (root / "fixture1/estimator.resume").exists()
    bare = root / "bare-state"
    run(2, fixture, "load", saved, bare, contains="VERIFIED_BYTES")
    assert not (bare / "estimator.resume").exists()
    run(4, estimator, *common, "--portable_load", bare, "--step_end", "1",
        contains="PORTABLE_RESTART circuit=unspecified step_start=0")
    # Sampler ignores damaged, mismatched or future estimator-only files.
    corrupt_resume = root / "corrupt-resume"; shutil.copytree(saved, corrupt_resume)
    (corrupt_resume / "estimator.resume").write_text("broken")
    output = run(2, sampler, "--portable_input", corrupt_resume, "--linefile", lines,
                 "--num_samples", "256", "--sampling_method", "1", "--bmps_bond_dim", "16")
    assert "estimator_next_step" not in output
    assert re.findall(r"^ sample \d+\.0: (\d+)$", output, re.M) == sampler_reference[2]
    run(2, estimator, *common, "--portable_load", corrupt_resume, "--step_end", "2",
        ok=False, contains="resume checksum mismatch")
    mismatch_resume = root / "mismatch-resume"; shutil.copytree(saved, mismatch_resume)
    change_resume(mismatch_resume, lambda s: re.sub(r'STATE "[a-f0-9]+"', 'STATE "wrong-state"', s))
    run(2, estimator, *common, "--portable_load", mismatch_resume, "--step_end", "2",
        ok=False, contains="resume state digest mismatch")
    # Width is derived from IDs, not the number of MPI ranks or a six-digit cap.
    wide = root / "wide-labels"
    run(2, fixture, "wide-save", wide)
    assert sorted(p.name for p in (wide / "vertices").iterdir()) == [
        "site-0000000.ten", "site-0000001.ten", "site-0000009.ten", "site-1000000.ten"]
    assert sorted(p.name for p in (wide / "messages").iterdir()) == [
        f"edge-{i:06d}-{direction}.ten" for i in range(3) for direction in ("forward", "reverse")]
    run(4, fixture, "load", wide, root / "wide-roundtrip", contains="VERIFIED_BYTES")
    if len(sys.argv) > 6:
        base = Path(sys.argv[6]); base.mkdir(parents=True, exist_ok=True)
        evidence = Path(tempfile.mkdtemp(prefix="accepted-", dir=base))
        for writers in (1, 2, 4): shutil.copytree(root / f"estimator{writers}", evidence / f"checkpoint-{writers}rank")
        for name in ("cluster.qasm", "inverse.qasm", "flip.qasm", "pauli.txt", "lines.txt"): shutil.copy2(root / name, evidence / name)
        for case, command, output in transcripts:
            (evidence / f"case-{case:02d}.log").write_text("COMMAND " + " ".join(command) + "\n" + output)
        print(f"EVIDENCE {evidence}")
    print(f"PASS {counter} portable checkpoint acceptance cases")
