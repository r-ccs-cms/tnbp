# kicked_ising_to_qasm2.py
from __future__ import annotations

import argparse
from typing import Iterable, Tuple, Set

from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import CouplingMap
from qiskit.qasm2 import dump as qasm2_dump, dumps as qasm2_dumps
from qiskit_ibm_runtime import QiskitRuntimeService

import math

def get_backend_and_cmap(backend_name: str):
    svc = QiskitRuntimeService()
    backend = svc.backend(backend_name)
    raw_cmap = backend.coupling_map or backend.configuration().coupling_map
    cmap = CouplingMap(couplinglist=raw_cmap)
    n_qubits = backend.num_qubits
    return backend, cmap, n_qubits


def undirected_edges_from_cmap(cmap: CouplingMap) -> Set[Tuple[int, int]]:
    edges: Set[Tuple[int, int]] = set()
    for (u, v) in cmap.get_edges():
        a, b = sorted((u, v))
        if a != b:
            edges.add((a, b))
    return edges

def greedy_matching_layers(edges: List[Tuple[int, int]]) -> List[List[Tuple[int, int]]]:
    """
    Decompose edges into layers where each layer is a matching (no shared vertices).
    Greedy construction; heavy-hex often fits in ~3 layers, but may vary.
    """
    remaining: Set[Tuple[int, int]] = set(edges)
    layers: List[List[Tuple[int, int]]] = []
    while remaining:
        used: Set[int] = set()
        layer: List[Tuple[int, int]] = []
        for (u, v) in sorted(remaining):
            if (u not in used) and (v not in used):
                layer.append((u, v))
                used.add(u)
                used.add(v)
        for e in layer:
            remaining.remove(e)
        layers.append(layer)
    return layers


def kicked_ising_layer(
    qc: QuantumCircuit,
    layers: List[List[Tuple[int, int]]],
    theta_x: float,
    theta_z: float,
        phi: float,
) -> int:

    n = qc.num_qubits
    count = 0

    # single-qubit kick
    for q in range(n):
        qc.rx(theta_x, q)
        count += 1
    for q in range(n):
        qc.rz(theta_z, q)
        count += 1

    # ZZ interation per parallelizable layer (cx-rz-cx)
    for layer in layers:
        for i, j in layer:
            qc.cx(i, j)
            qc.rz(phi, j)
            qc.cx(i, j)
            count += 3
        # barrier between parallel layers
        qc.barrier(*qc.qubits)

    return count


def build_kicked_ising_circuit(
    n_qubits: int,
    edges: Iterable[Tuple[int, int]],
    steps: int,
    theta_x: float,
    theta_z: float,
        phi: float,
    measure: bool = False,
):
    qc = QuantumCircuit(n_qubits, n_qubits if measure else 0)
    step_counts = []

        # compute parallelizable layers from edges
    layers = greedy_matching_layers(edges)
    for step in range(steps):
        c = kicked_ising_layer(qc, layers, theta_x, theta_z, phi)
        step_counts.append(c)

    if measure:
        qc.barrier()
        qc.measure(range(n_qubits), range(n_qubits))

    return qc, step_counts


def main():
    ap = argparse.ArgumentParser(description="Dump QASM file for kicked Ising on backend topology")
    ap.add_argument("--backend", type=str, default="ibm_kobe", help="IBM backend name (e.g., 'ibm_kobe').")
    ap.add_argument("--steps", type=int, default=1, help="Number of Floquet steps")
    ap.add_argument("--hx", type=float, default=0.5, help="Coefficient for sum_i X_i (default 0.25)")
    ap.add_argument("--hz", type=float, default=0.0, help="Coefficient for sum_i Z_i (default 0.0)")
    ap.add_argument("--jz", type=float, default=0.5, help="Coefficient for sum_(i,j in edges) Z_i Z_j (default 0.25)")
    ap.add_argument("--output", type=str, default="kicked_ising.qasm", help="Output path (qasm)")
    args = ap.parse_args()
    
    backend_name = args.backend
    steps = args.steps
    theta_x = args.hx * math.pi
    theta_z = args.hz * math.pi
    phi = args.jz * math.pi
    add_measure = False
    qasm_out = args.output

    backend, cmap, n_qubits = get_backend_and_cmap(backend_name)
    bonds = undirected_edges_from_cmap(cmap)

    qc, step_counts = build_kicked_ising_circuit(
        n_qubits=n_qubits,
        edges=sorted(bonds),
        steps=steps,
        theta_x=theta_x,
        theta_z=theta_z,
        phi=phi,
        measure=add_measure,
    )

    basis = ["u", "cx"]
    routed = transpile(qc, coupling_map=cmap, basis_gates=basis, optimization_level=1)

    # try:
    # Qiskit 1.x 系の推奨：プリセット PM を使う（バリアは尊重される）
    #    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    #    pm = generate_preset_pass_manager(
    #        optimization_level=1,           # 0〜1 を推奨（順序保全重視）
    #        target=backend.target if backend else None,
    #        basis_gates=basis,
    #        layout_method="sabre",
    #        routing_method="sabre",
    #        # 注意：RemoveBarriers を入れないプリセット（デフォルトでもバリアは尊重されます）
    #    )
    #    routed = pm.run(qc)
    #except Exception:
    # 旧 API 互換ルート
    #    routed = transpile(
    #        qc,
    #        backend=backend if backend else None,
    #        coupling_map=cmap if cmap else None,
    #        basis_gates=basis,
    #        layout_method="sabre",
    #        routing_method="sabre",
    #        optimization_level=1,  # 0 でも OK
    #        seed_transpiler=42,
    #    )    

    # QASM 出力
    qasm2_dump(routed, qasm_out)

    # === 追加: step ごとのゲート数表示 ===
    print("=== Gate counts per step (before transpile) ===")
    for i, c in enumerate(step_counts, 1):
        print(f"Step {i}: {c} gates")

    print(f"\n[OK] Wrote OpenQASM 2.0 to: {qasm_out}")


if __name__ == "__main__":
    main()
