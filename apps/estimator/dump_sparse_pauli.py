#!/usr/bin/env python3
# Build (transverse+longitudinal-field Ising) on a backend topology and dump as "real imag PauliString"

import argparse
from typing import List, Tuple, Optional
from collections import defaultdict

from qiskit.quantum_info import SparsePauliOp

# Optional (only needed when --backend is IBM device)
# If you don't use IBM backends, you can remove these imports.
try:
    from qiskit_ibm_runtime import QiskitRuntimeService
except Exception:
    QiskitRuntimeService = None  # type: ignore


def one_site_label(n: int, i: int, axis: str) -> str:
    # Qiskit: rightmost char is qubit 0
    s = ['I'] * n
    s[n - 1 - i] = axis
    return ''.join(s)


def two_site_label(n: int, i: int, j: int, axis_i: str, axis_j: str) -> str:
    s = ['I'] * n
    s[n - 1 - i] = axis_i
    s[n - 1 - j] = axis_j
    return ''.join(s)


def unique_undirected_edges(edges: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    # normalize (min,max) and deduplicate
    seen = set()
    out = []
    for a, b in edges:
        if a == b:
            continue
        e = (a, b) if a < b else (b, a)
        if e not in seen:
            seen.add(e)
            out.append(e)
    return out


def get_backend_and_topology(
    backend_name: Optional[str],
    ibm_instance: Optional[str],
) -> Tuple[Optional[object], List[Tuple[int, int]], int]:
    """
    Returns: (backend_or_None, edges, n_qubits)
    - If backend_name is None: uses a fully-connected mock topology (requires --n-qubits)
    - If IBM backend name is given: queries the service and extracts coupling edges
    """
    # No backend given → caller must supply --n-qubits and (optionally) --edge-file
    if backend_name is None:
        return None, [], -1

    if QiskitRuntimeService is None:
        raise RuntimeError(
            "qiskit_ibm_runtime is not available but --backend was specified. "
            "Install it or run without --backend."
        )

    service = QiskitRuntimeService(instance=ibm_instance) if ibm_instance else QiskitRuntimeService()
    backend = service.backend(backend_name)

    # Try multiple ways to get coupling map
    edges = None

    # 1) direct attribute
    try:
        cm = getattr(backend, "coupling_map", None)
        if cm:
            edges = [(int(a), int(b)) for a, b in cm]
    except Exception:
        pass

    # 2) configuration().coupling_map
    if edges is None:
        try:
            conf = backend.configuration()
            cm = getattr(conf, "coupling_map", None)
            if cm:
                edges = [(int(a), int(b)) for a, b in cm]
        except Exception:
            pass

    # 3) target.build_coupling_map()
    if edges is None:
        try:
            target = backend.target
            cm = target.build_coupling_map()
            if cm:
                edges = [(int(a), int(b)) for a, b in cm]
        except Exception:
            pass

    if edges is None:
        raise RuntimeError(f"Could not obtain coupling map from backend '{backend_name}'")

    edges = unique_undirected_edges(edges)
    # n_qubits
    try:
        n = int(getattr(backend, "num_qubits"))
    except Exception:
        try:
            n = int(backend.configuration().num_qubits)
        except Exception:
            # Fall back to max index in edges + 1
            n = max(max(a, b) for a, b in edges) + 1

    return backend, edges, n


def build_ising_on_topology(n: int, edges: List[Tuple[int, int]], hx: float, hz: float, jz: float) -> SparsePauliOp:
    """
    H = hx * sum_i X_i + hz * sum_i Z_i + jz * sum_(i,j in edges) Z_i Z_j
    """
    terms = []

    if hx != 0.0:
        for i in range(n):
            terms.append((one_site_label(n, i, 'X'), float(hx)))

    if hz != 0.0:
        for i in range(n):
            terms.append((one_site_label(n, i, 'Z'), float(hz)))

    if jz != 0.0:
        for (i, j) in edges:
            terms.append((two_site_label(n, i, j, 'Z', 'Z'), float(jz)))

    if not terms:
        raise ValueError("No terms were added. Set at least one of hx, hz, jz non-zero.")
    return SparsePauliOp.from_list(terms)


def dump_txt(op, path):
    """
    行形式: 'real imag PauliString'
    - 恒等 (I...I) はスキップ
    - 同じ文字列は係数を合算
    """
    acc = defaultdict(complex)
    for p, c in zip(op.paulis, op.coeffs):
        s = p.to_label()          # ← これがポイント！ str(p) はダメ
        if all(ch == 'I' for ch in s):
            continue
        acc[s] += complex(c)

    EPS = 1e-12
    with open(path, "w") as f:
        for s, c in acc.items():
            if abs(c.real) < EPS and abs(c.imag) < EPS:
                continue
            f.write(f"{c.real} {c.imag} {s}\n")


def main():
    ap = argparse.ArgumentParser(description="Dump (transverse+longitudinal-field Ising) on backend topology")
    ap.add_argument("--backend", type=str, default=None, help="IBM backend name (e.g., 'ibm_osaka'). If omitted, you must provide --n-qubits and --edge-file or rely on fully-connected if --fully-connected is set.")
    ap.add_argument("--ibm-instance", type=str, default=None, help="IBM Cloud instance string if needed (e.g., 'hub/group/project').")
    ap.add_argument("--n-qubits", type=int, default=None, help="Number of qubits (required if --backend not given).")
    ap.add_argument("--edge-file", type=str, default=None, help="Optional text file of edges: each line 'i j'.")
    ap.add_argument("--fully-connected", action="store_true", help="If set (and --backend not given), build fully-connected edges over n qubits.")
    ap.add_argument("--hx", type=float, default=1.0, help="Coefficient for sum_i X_i (default 1.0)")
    ap.add_argument("--hz", type=float, default=1.0, help="Coefficient for sum_i Z_i (default 1.0)")
    ap.add_argument("--jz", type=float, default=1.0, help="Coefficient for sum_(i,j in edges) Z_i Z_j (default 1.0)")
    ap.add_argument("--out", type=str, default="ising_topology.txt", help="Output path (txt)")
    args = ap.parse_args()

    # Determine topology & n
    if args.backend:
        _, edges, n = get_backend_and_topology(args.backend, args.ibm_instance)
    else:
        if args.n_qubits is None:
            raise ValueError("When --backend is not provided, you must set --n-qubits.")
        n = args.n_qubits
        edges = []
        if args.edge_file:
            with open(args.edge_file, "r") as ef:
                for line in ef:
                    line = line.strip()
                    if not line or line.startswith("#") or line.startswith("//"):
                        continue
                    a_str, b_str = line.split()
                    edges.append((int(a_str), int(b_str)))
            edges = unique_undirected_edges(edges)
        elif args.fully_connected:
            edges = unique_undirected_edges([(i, j) for i in range(n) for j in range(i + 1, n)])
        else:
            raise ValueError("No topology provided. Use --edge-file or --fully-connected (or specify --backend).")

    op = build_ising_on_topology(n, edges, args.hx, args.hz, args.jz)
    dump_txt(op, args.out)
    print(f"Wrote {len(op.paulis)} terms to {args.out} (n_qubits={n}, edges={len(edges)})")


if __name__ == "__main__":
    main()
