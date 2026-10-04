#!/usr/bin/env python3

import sys
import math
import random
import argparse
import matplotlib.pyplot as plt


def read_edge_list(filename):
    edges = []
    vertices = set()
    n_sites = None

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("#"):
                # optional: "# n_sites = 10" を拾う
                s = line[1:].strip()
                if "=" in s:
                    lhs, rhs = s.split("=", 1)
                    lhs = lhs.strip()
                    rhs = rhs.strip()
                    if lhs == "n_sites":
                        try:
                            n_sites = int(rhs)
                        except ValueError:
                            pass
                continue

            parts = line.split()
            if len(parts) != 2:
                raise ValueError(f"invalid line: {line}")

            u = int(parts[0])
            v = int(parts[1])

            edges.append((u, v))
            vertices.add(u)
            vertices.add(v)

    if n_sites is not None:
        for i in range(n_sites):
            vertices.add(i)

    vertices = sorted(vertices)
    return vertices, edges


def initial_positions(vertices, seed=12345):
    rng = random.Random(seed)
    pos = {}
    n = len(vertices)

    if n == 0:
        return pos

    # 円周上 + 少しノイズ
    for k, v in enumerate(vertices):
        theta = 2.0 * math.pi * k / n
        r = 1.0 + 0.05 * rng.random()
        x = r * math.cos(theta)
        y = r * math.sin(theta)
        pos[v] = [x, y]

    return pos


def spring_layout(vertices, edges, seed=12345, iterations=200, k=None):
    """
    簡易 force-directed layout
    networkx なしで動くように自前実装
    """
    pos = initial_positions(vertices, seed=seed)
    n = len(vertices)

    if n == 0:
        return pos

    index_of = {v: i for i, v in enumerate(vertices)}

    if k is None:
        k = 1.0 / math.sqrt(max(n, 1))

    temperature = 0.1

    for _ in range(iterations):
        disp = {v: [0.0, 0.0] for v in vertices}

        # 反発力
        for i in range(n):
            vi = vertices[i]
            xi, yi = pos[vi]
            for j in range(i + 1, n):
                vj = vertices[j]
                xj, yj = pos[vj]

                dx = xi - xj
                dy = yi - yj
                dist2 = dx * dx + dy * dy + 1.0e-9
                dist = math.sqrt(dist2)

                force = (k * k) / dist

                fx = force * dx / dist
                fy = force * dy / dist

                disp[vi][0] += fx
                disp[vi][1] += fy
                disp[vj][0] -= fx
                disp[vj][1] -= fy

        # 引力
        for u, v in edges:
            xu, yu = pos[u]
            xv, yv = pos[v]

            dx = xu - xv
            dy = yu - yv
            dist2 = dx * dx + dy * dy + 1.0e-9
            dist = math.sqrt(dist2)

            force = dist2 / k

            fx = force * dx / dist
            fy = force * dy / dist

            disp[u][0] -= fx
            disp[u][1] -= fy
            disp[v][0] += fx
            disp[v][1] += fy

        # 更新
        for v in vertices:
            dx, dy = disp[v]
            d = math.sqrt(dx * dx + dy * dy)
            if d > 0.0:
                scale = min(d, temperature) / d
                pos[v][0] += dx * scale
                pos[v][1] += dy * scale

        temperature *= 0.98

    return pos


def draw_graph(vertices, edges, pos, output=None, title=None,
               node_size=0.08, font_size=10, edge_width=1.0):
    fig, ax = plt.subplots(figsize=(8, 8))

    # edges
    for u, v in edges:
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        ax.plot([x1, x2], [y1, y2],
                linewidth=edge_width,
                color="#000000",
                zorder=1)

    # nodes
    for v in vertices:
        x, y = pos[v]
        circle = plt.Circle((x, y),
                            node_size,
                            fill=True,
                            facecolor="#eeffee",
                            edgecolor="#000000",
                            linewidth=edge_width,
                            zorder=2)
        ax.add_patch(circle)

        ax.text(
            x, y, str(v),
            ha="center", va="center",
            fontsize=font_size,
            color="#000000",
            zorder=3
        )

    if title is not None:
        ax.set_title(title)

    ax.set_aspect("equal")
    ax.axis("off")

    # 余白を少し取る
    if vertices:
        xs = [pos[v][0] for v in vertices]
        ys = [pos[v][1] for v in vertices]
        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)
        margin = 0.2
        ax.set_xlim(xmin - margin, xmax + margin)
        ax.set_ylim(ymin - margin, ymax + margin)

    plt.tight_layout()

    if output is not None:
        plt.savefig(output, dpi=200)
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description="Draw graph from edge list.")
    parser.add_argument("input", help="input edge-list file")
    parser.add_argument("-o", "--output", default=None, help="output image file")
    parser.add_argument("--seed", type=int, default=12345, help="random seed")
    parser.add_argument("--iter", type=int, default=200, help="layout iterations")
    parser.add_argument("--title", default=None, help="plot title")
    parser.add_argument("--node-size", type=float, default=0.08, help="node radius")
    parser.add_argument("--font-size", type=int, default=10, help="label font size")
    parser.add_argument("--edge-width", type=float, default=1.0,
                        help="edge line width")
    parser.add_argument("--auto-scale", action="store_true",
                        help="auto adjust sizes based on number of nodes")
    args = parser.parse_args()

    vertices, edges = read_edge_list(args.input)
    pos = spring_layout(vertices, edges, seed=args.seed, iterations=args.iter)

    title = args.title if args.title is not None else args.input

    n = len(vertices)

    node_size = args.node_size
    font_size = args.font_size
    edge_width = args.edge_width

    if args.auto_scale and n > 0:
        scale = 1.0 / (n ** 0.5)
        node_size *= 2.0 * scale
        edge_width *= 1.5 * scale
        font_size = max(4, int(font_size * scale * 2.5))
    
    draw_graph(
        vertices, edges, pos,
        output=args.output,
        title=title,
        node_size=node_size,
        font_size=font_size,
        edge_width=edge_width
    )


if __name__ == "__main__":
    main()
