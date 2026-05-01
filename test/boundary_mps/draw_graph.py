#!/usr/bin/env python3
"""
Draw an edge-list graph with a tunable force-directed layout.

The script is still for general graphs, but adds options that are useful for
sequences of gradually changing graphs:
  * deterministic non-random initial positions: --init circle/grid/path/random
  * read/write positions: --pos-in previous.pos --pos-out current.pos
  * force tuning: --k, --edge-length, --repulsion, --attraction, --gravity
  * optional multiple random restarts for a single static graph: --restarts N

Input format:
  lines with "u v" edges; comment lines beginning with # are ignored.
  If a comment contains "# n_sites = N", isolated vertices 0..N-1 are included.
"""

import argparse
import math
import random
from pathlib import Path



def read_edge_list(filename):
    edges = []
    vertices = set()
    n_sites = None

    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                s = line[1:].strip()
                if "=" in s:
                    lhs, rhs = s.split("=", 1)
                    if lhs.strip() == "n_sites":
                        try:
                            n_sites = int(rhs.strip())
                        except ValueError:
                            pass
                continue
            parts = line.split()
            if len(parts) != 2:
                raise ValueError(f"invalid line: {line}")
            u, v = int(parts[0]), int(parts[1])
            edges.append((u, v))
            vertices.add(u)
            vertices.add(v)

    if n_sites is not None:
        vertices.update(range(n_sites))

    return sorted(vertices), edges


def connected_components(vertices, edges):
    adj = {v: [] for v in vertices}
    for u, v in edges:
        if u in adj and v in adj:
            adj[u].append(v)
            adj[v].append(u)
    seen = set()
    comps = []
    for s in vertices:
        if s in seen:
            continue
        stack = [s]
        seen.add(s)
        comp = []
        while stack:
            v = stack.pop()
            comp.append(v)
            for w in adj[v]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        comps.append(sorted(comp))
    comps.sort(key=lambda c: (-len(c), c[0] if c else 0))
    return comps


def initial_positions(vertices, seed=12345, mode="circle", jitter=0.02):
    rng = random.Random(seed)
    n = len(vertices)
    pos = {}
    if n == 0:
        return pos

    if mode == "random":
        for v in vertices:
            pos[v] = [rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0)]

    elif mode == "grid":
        # Generic, label-order based initializer.  This is not honeycomb-specific,
        # but works well when vertex labels roughly follow the geometry.
        nx = math.ceil(math.sqrt(n))
        for k, v in enumerate(vertices):
            x = k % nx
            y = -(k // nx)
            pos[v] = [float(x), float(y)]

    elif mode == "path":
        # Good for chain-like or label-ordered graphs.
        for k, v in enumerate(vertices):
            pos[v] = [float(k), 0.0]

    else:  # circle
        for k, v in enumerate(vertices):
            theta = 2.0 * math.pi * k / n
            pos[v] = [math.cos(theta), math.sin(theta)]

    if jitter > 0:
        for v in vertices:
            pos[v][0] += jitter * (rng.random() - 0.5)
            pos[v][1] += jitter * (rng.random() - 0.5)
    return pos


def read_positions(filename, vertices, fallback_pos):
    pos = {v: list(fallback_pos[v]) for v in vertices}
    p = Path(filename)
    if not p.exists():
        raise FileNotFoundError(filename)
    with open(p, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            v = int(parts[0])
            if v in pos:
                pos[v] = [float(parts[1]), float(parts[2])]
    return pos


def write_positions(filename, vertices, pos):
    with open(filename, "w", encoding="utf-8") as f:
        f.write("# v x y\n")
        for v in vertices:
            x, y = pos[v]
            f.write(f"{v} {x:.17g} {y:.17g}\n")


def normalize_positions(vertices, pos, scale=1.0):
    if not vertices:
        return pos
    xs = [pos[v][0] for v in vertices]
    ys = [pos[v][1] for v in vertices]
    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    maxr = 0.0
    for v in vertices:
        pos[v][0] -= cx
        pos[v][1] -= cy
        maxr = max(maxr, math.hypot(pos[v][0], pos[v][1]))
    if maxr > 0:
        s = scale / maxr
        for v in vertices:
            pos[v][0] *= s
            pos[v][1] *= s
    return pos


def layout_energy(vertices, edges, pos, edge_length=1.0):
    e = 0.0
    for u, v in edges:
        dx = pos[u][0] - pos[v][0]
        dy = pos[u][1] - pos[v][1]
        d = math.hypot(dx, dy)
        e += (d - edge_length) ** 2
    return e


def spring_layout(vertices, edges, seed=12345, iterations=800, init="circle", pos0=None,
                  k=None, edge_length=1.0, repulsion=1.0, attraction=1.0,
                  gravity=0.03, temperature=0.25, cooling=0.995):
    n = len(vertices)
    if n == 0:
        return {}

    pos = {v: list(pos0[v]) for v in vertices} if pos0 is not None else initial_positions(vertices, seed, init)
    pos = normalize_positions(vertices, pos, scale=max(1.0, math.sqrt(n)))

    if k is None:
        # Larger than the classical 1/sqrt(n) default.  This avoids collapsing
        # sparse lattices into a hairball.
        k = edge_length * math.sqrt(max(n, 1)) / 2.0

    temp = temperature * max(1.0, math.sqrt(n))

    for _ in range(iterations):
        disp = {v: [0.0, 0.0] for v in vertices}

        # Pair repulsion.
        for i, vi in enumerate(vertices):
            xi, yi = pos[vi]
            for vj in vertices[i + 1:]:
                xj, yj = pos[vj]
                dx = xi - xj
                dy = yi - yj
                dist2 = dx * dx + dy * dy + 1.0e-12
                dist = math.sqrt(dist2)
                force = repulsion * k * k / dist
                fx = force * dx / dist
                fy = force * dy / dist
                disp[vi][0] += fx
                disp[vi][1] += fy
                disp[vj][0] -= fx
                disp[vj][1] -= fy

        # Edge attraction with a nonzero preferred edge length.
        for u, v in edges:
            xu, yu = pos[u]
            xv, yv = pos[v]
            dx = xu - xv
            dy = yu - yv
            dist = math.hypot(dx, dy) + 1.0e-12
            # Hooke-like force; zero around edge_length.
            force = attraction * (dist - edge_length)
            fx = force * dx / dist
            fy = force * dy / dist
            disp[u][0] -= fx
            disp[u][1] -= fy
            disp[v][0] += fx
            disp[v][1] += fy

        # Weak gravity to keep components on the page without crushing them.
        if gravity != 0.0:
            for v in vertices:
                disp[v][0] -= gravity * pos[v][0]
                disp[v][1] -= gravity * pos[v][1]

        for v in vertices:
            dx, dy = disp[v]
            d = math.hypot(dx, dy)
            if d > 0.0:
                step = min(d, temp) / d
                pos[v][0] += dx * step
                pos[v][1] += dy * step

        temp *= cooling

    return normalize_positions(vertices, pos, scale=max(1.0, math.sqrt(n)))


def draw_graph(vertices, edges, pos, output=None, title=None, node_size=0.08,
               font_size=10, edge_width=1.0, labels=True, figsize="8x8",
               margin=0.3):
    import matplotlib.pyplot as plt
    if "x" in figsize:
        w, h = [float(x) for x in figsize.lower().split("x", 1)]
    else:
        w = h = float(figsize)
    fig, ax = plt.subplots(figsize=(w, h))

    for u, v in edges:
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        ax.plot([x1, x2], [y1, y2], linewidth=edge_width, color="black", zorder=1)

    for v in vertices:
        x, y = pos[v]
        circ = plt.Circle((x, y), node_size, fill=True, facecolor="#eeffee",
                          edgecolor="black", linewidth=edge_width, zorder=2)
        ax.add_patch(circ)
        if labels:
            ax.text(x, y, str(v), ha="center", va="center", fontsize=font_size,
                    color="black", zorder=3)

    if title:
        ax.set_title(title)
    ax.set_aspect("equal")
    ax.axis("off")

    if vertices:
        xs = [pos[v][0] for v in vertices]
        ys = [pos[v][1] for v in vertices]
        ax.set_xlim(min(xs) - margin, max(xs) + margin)
        ax.set_ylim(min(ys) - margin, max(ys) + margin)

    plt.tight_layout()
    if output:
        plt.savefig(output, dpi=200)
    else:
        plt.show()
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Draw a general graph from an edge list.")
    ap.add_argument("input")
    ap.add_argument("-o", "--output", default=None)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--iter", type=int, default=800)
    ap.add_argument("--init", choices=["circle", "grid", "path", "random"], default="circle")
    ap.add_argument("--pos-in", default=None, help="read initial positions from file")
    ap.add_argument("--pos-out", default=None, help="write final positions to file")
    ap.add_argument("--restarts", type=int, default=1, help="try several seeds and keep the lowest edge-length energy")
    ap.add_argument("--k", type=float, default=None)
    ap.add_argument("--edge-length", type=float, default=1.0)
    ap.add_argument("--repulsion", type=float, default=1.0)
    ap.add_argument("--attraction", type=float, default=1.0)
    ap.add_argument("--gravity", type=float, default=0.03)
    ap.add_argument("--temperature", type=float, default=0.25)
    ap.add_argument("--cooling", type=float, default=0.995)
    ap.add_argument("--title", default=None)
    ap.add_argument("--node-size", type=float, default=0.08)
    ap.add_argument("--font-size", type=int, default=10)
    ap.add_argument("--edge-width", type=float, default=1.0)
    ap.add_argument("--figsize", default="8x8")
    ap.add_argument("--margin", type=float, default=0.3)
    ap.add_argument("--no-labels", action="store_true")
    ap.add_argument("--auto-scale", action="store_true")
    args = ap.parse_args()

    vertices, edges = read_edge_list(args.input)
    fallback = initial_positions(vertices, args.seed, args.init)
    pos0 = read_positions(args.pos_in, vertices, fallback) if args.pos_in else fallback

    best_pos = None
    best_e = None
    restarts = max(1, args.restarts)
    for r in range(restarts):
        this_pos0 = pos0 if args.pos_in else initial_positions(vertices, args.seed + r, args.init)
        pos = spring_layout(vertices, edges, seed=args.seed + r, iterations=args.iter,
                            init=args.init, pos0=this_pos0, k=args.k,
                            edge_length=args.edge_length, repulsion=args.repulsion,
                            attraction=args.attraction, gravity=args.gravity,
                            temperature=args.temperature, cooling=args.cooling)
        e = layout_energy(vertices, edges, pos, args.edge_length)
        if best_e is None or e < best_e:
            best_e = e
            best_pos = pos

    if args.pos_out:
        write_positions(args.pos_out, vertices, best_pos)

    node_size = args.node_size
    font_size = args.font_size
    edge_width = args.edge_width
    if args.auto_scale and vertices:
        s = 1.0 / math.sqrt(len(vertices))
        node_size *= 2.2 * s
        edge_width *= 1.8 * s
        font_size = max(4, int(font_size * 2.8 * s))

    draw_graph(vertices, edges, best_pos, output=args.output,
               title=args.title if args.title is not None else args.input,
               node_size=node_size, font_size=font_size, edge_width=edge_width,
               labels=not args.no_labels, figsize=args.figsize, margin=args.margin)


if __name__ == "__main__":
    main()
