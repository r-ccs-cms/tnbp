#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from typing import List, Tuple, Set


def normalize_edges(coupling_map: List[List[int]]) -> List[Tuple[int, int]]:
    """Normalize (u, v) as u < v, drop self-loops/dupes, and sort."""
    s: Set[Tuple[int, int]] = set()
    for e in coupling_map:
        if len(e) != 2:
            continue
        u, v = int(e[0]), int(e[1])
        if u == v:
            continue
        a, b = (u, v) if u < v else (v, u)
        s.add((a, b))
    return sorted(s)


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


def to_cpp_pairs_list(edges: List[Tuple[int, int]]) -> str:
    """
    Make a C++ initializer list string: { {u,v}, {x,y}, ... } with line wrapping.
    """
    items = [f"{{{u}, {v}}}" for (u, v) in edges]
    lines = []
    for i in range(0, len(items), 10):
        lines.append(", ".join(items[i : i + 10]))
    body = (",\n        ".join(lines)) if lines else ""
    return "{ " + ("\n        " + body + "\n    " if body else "") + "}"


def make_safe_identifier(name: str) -> str:
    """
    Convert arbitrary backend name into a safe C++ identifier (letters, digits, '_').
    e.g., 'ibm-marrakesh' -> 'ibm_marrakesh'
    """
    return re.sub(r"[^0-9A-Za-z_]", "_", name)


def emit_cpp_backend_functions(
    backend_name_raw: str,
    n_qubits: int,
    edges: List[Tuple[int, int]],
    layers: List[List[Tuple[int, int]]],
) -> str:
    """
    Emit a full C++ header:
      - file header comment
      - include guard (#ifndef / #define / #endif)
      - namespace tnbp { ... }
      - two functions:
          1) bond_<backend>()
          2) parallel_bond_<backend>()
    """

    safe_name = make_safe_identifier(backend_name_raw)  # e.g., ibm_marrakesh
    guard_name = f"TNBP_DEVICE_{safe_name.upper()}_H"

    header = f"""//// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/device/{safe_name}.h
@brief coupling map for {backend_name_raw}
*/
#ifndef {guard_name}
#define {guard_name}

#include <vector>
#include <utility>

namespace tnbp {{

"""

    fn1 = f"""/// Coupling map for {backend_name_raw} (n_qubits = {n_qubits})
/// Edges are undirected and normalized as (u<v).
inline std::vector<std::pair<int,int>> bond_{safe_name}() {{
    return std::vector<std::pair<int,int>>{to_cpp_pairs_list(edges)};
}}
"""

    layer_bodies = []
    for layer in layers:
        layer_bodies.append("std::vector<std::pair<int,int>>" + to_cpp_pairs_list(layer))
    layers_joined = ",\n        ".join(layer_bodies) if layer_bodies else ""

    fn2 = f"""
/// Parallelizable bond layers (matchings) for {backend_name_raw}.
/// Edges in the same inner vector can be executed in parallel without qubit conflicts.
/// The decomposition is greedy; the number of layers may vary by topology.
inline std::vector<std::vector<std::pair<int,int>>> parallel_bond_{safe_name}() {{
    return std::vector<std::vector<std::pair<int,int>>> {{
        {layers_joined}
    }};
}}
"""

    footer = f"""
}} // namespace tnbp

#endif // {guard_name}
"""
    return header + fn1 + fn2 + footer


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Emit C++ header that hard-codes a backend's coupling map and "
            "its parallelizable layers (matchings). Prints to stdout."
        )
    )
    parser.add_argument(
        "--backend",
        default="ibm_marrakesh",
        help="Backend name (default: ibm_marrakesh)",
    )
    parser.add_argument(
        "--offline-json",
        default=None,
        help=(
            "Optional path to JSON with fields "
            "'n_qubits' (int) and 'coupling_map' ([[u,v],...]). "
            "Use this when IBM cloud access is unavailable."
        ),
    )
    args = parser.parse_args()

    if args.offline_json:
        import json
        import pathlib

        data = json.loads(pathlib.Path(args.offline_json).read_text())
        coupling_map = data["coupling_map"]
        n_qubits = int(data.get("n_qubits", 0))
        backend_name_raw = args.backend
    else:
        # Online query via qiskit-ibm-runtime (requires prior token setup)
        from qiskit_ibm_runtime import QiskitRuntimeService

        svc = QiskitRuntimeService()
        backend = svc.backend(args.backend)
        cfg = backend.configuration()
        coupling_map = cfg.coupling_map
        n_qubits = int(cfg.n_qubits)
        backend_name_raw = args.backend

    edges = normalize_edges(coupling_map)
    layers = greedy_matching_layers(edges)

    cpp = emit_cpp_backend_functions(backend_name_raw, n_qubits, edges, layers)
    print(cpp)


if __name__ == "__main__":
    main()
