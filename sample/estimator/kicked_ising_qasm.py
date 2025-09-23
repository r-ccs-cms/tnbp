# kicked_ising_to_qasm2.py
from __future__ import annotations
from typing import Iterable, Tuple, Set

from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import CouplingMap
from qiskit.qasm2 import dump as qasm2_dump, dumps as qasm2_dumps
from qiskit_ibm_runtime import QiskitRuntimeService


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


def kicked_ising_layer(
    qc: QuantumCircuit,
    edges: Iterable[Tuple[int, int]],
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

    # ZZ interation (cx-rz-cx)
    
    for i, j in edges:
        qc.cx(i, j)
        qc.rz(phi, j)
        qc.cx(i, j)
        count += 3

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

    for step in range(steps):
        c = kicked_ising_layer(qc, edges, theta_x, theta_z, phi)
        step_counts.append(c)

    if measure:
        qc.barrier()
        qc.measure(range(n_qubits), range(n_qubits))

    return qc, step_counts


def main():
    backend_name = "ibm_kobe"
    steps = 1
    theta_x = 0.9 * 3.14159265
    theta_z = 0.0
    phi = 1.0 * 0.5 * 3.14159265
    add_measure = False
    qasm_out = "kicked_ising.qasm"

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

    # QASM 出力
    qasm2_dump(routed, qasm_out)

    # === 追加: step ごとのゲート数表示 ===
    print("=== Gate counts per step (before transpile) ===")
    for i, c in enumerate(step_counts, 1):
        print(f"Step {i}: {c} gates")

    print(f"\n[OK] Wrote OpenQASM 2.0 to: {qasm_out}")


if __name__ == "__main__":
    main()
