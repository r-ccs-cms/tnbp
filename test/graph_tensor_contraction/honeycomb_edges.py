#!/usr/bin/env python3

import argparse


def site_index(ia, ib, sub, La, Lb):
    """
    sub = 0 (A), 1 (B)

    ラベルは a方向を最速に進める:
        (ib fixed, ia = 0..La-1) を先に回す
    各ユニットセル内では A, B の順。
    """
    return 2 * (ib * La + ia) + sub


def make_edge(i, j):
    return (i, j) if i < j else (j, i)


def honeycomb_edges_open(La, Lb):
    """
    Open boundary condition の honeycomb lattice の辺リストを返す。
    各 A サイトから 3 本（境界では減る）の結合を張る:
        A(ia,ib) -- B(ia,ib)
        A(ia,ib) -- B(ia-1,ib)   if ia > 0
        A(ia,ib) -- B(ia,ib-1)   if ib > 0
    """
    if La <= 0 or Lb <= 0:
        raise ValueError("La and Lb must be positive")

    edges = []

    for ib in range(Lb):
        for ia in range(La):
            a = site_index(ia, ib, 0, La, Lb)  # A(ia,ib)

            # 1) same unit cell
            b0 = site_index(ia, ib, 1, La, Lb)
            edges.append(make_edge(a, b0))

            # 2) neighboring cell in -a direction
            if ia > 0:
                b1 = site_index(ia - 1, ib, 1, La, Lb)
                edges.append(make_edge(a, b1))

            # 3) neighboring cell in -b direction
            if ib > 0:
                b2 = site_index(ia, ib - 1, 1, La, Lb)
                edges.append(make_edge(a, b2))

    edges.sort()
    return edges


def main():
    parser = argparse.ArgumentParser(
        description="Generate open-boundary honeycomb lattice edge list."
    )
    parser.add_argument("La", type=int, help="number of unit cells along a")
    parser.add_argument("Lb", type=int, help="number of unit cells along b")
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="output filename (default: stdout)"
    )
    args = parser.parse_args()

    La = args.La
    Lb = args.Lb

    edges = honeycomb_edges_open(La, Lb)
    n_sites = 2 * La * Lb

    out = open(args.output, "w") if args.output is not None else None
    f = out if out is not None else __import__("sys").stdout

    try:
        print(f"# honeycomb lattice, open boundary", file=f)
        print(f"# La = {La}", file=f)
        print(f"# Lb = {Lb}", file=f)
        print(f"# n_sites = {n_sites}", file=f)
        print(f"# n_edges = {len(edges)}", file=f)
        print(f"# format: u v", file=f)

        for u, v in edges:
            print(f"{u} {v}", file=f)
    finally:
        if out is not None:
            out.close()


if __name__ == "__main__":
    main()