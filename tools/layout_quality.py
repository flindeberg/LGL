#!/usr/bin/env python3
# Layout quality metrics for comparing lglayout versions (developer tool, not
# part of the pipeline). Usage: layout_quality.py <graph.lgl> <coords>...
# One line per coords file: nodes, NaN/Inf count, extent, median/p90 edge
# length, coefficient of variation of edge length, share of nodes with a
# neighbour closer than 0.02 / 0.1, mean nearest-neighbour distance (< 0.5).
import sys, math, statistics
from collections import defaultdict

def read_lgl(path):
    edges = []
    cur = None
    with open(path) as f:
        for line in f:
            p = line.split()
            if not p:
                continue
            if p[0] == '#':
                cur = p[1]
            else:
                edges.append((cur, p[0]))
    return edges

def pct(xs, q):
    return xs[min(len(xs) - 1, int(q * len(xs)))]

edges = read_lgl(sys.argv[1])
print("file nodes bad extent_x extent_y med_edge p90_edge cv_edge overlap@0.02 overlap@0.1 mean_nn")
for c in sys.argv[2:]:
    pos = {}
    bad = 0
    with open(c) as f:
        for line in f:
            p = line.split()
            x, y = float(p[1]), float(p[2])
            if not (math.isfinite(x) and math.isfinite(y)):
                bad += 1
            pos[p[0]] = (x, y)
    xs = [v[0] for v in pos.values()]; ys = [v[1] for v in pos.values()]
    L = sorted(math.dist(pos[a], pos[b]) for a, b in edges)
    mean = statistics.fmean(L); sd = statistics.pstdev(L)
    # nearest-neighbour stats with a grid hash
    cell = 0.5
    grid = defaultdict(list)
    for k, (x, y) in pos.items():
        grid[(int(math.floor(x / cell)), int(math.floor(y / cell)))].append((x, y))
    o1 = o2 = 0; nn_sum = 0.0; nn_cnt = 0
    for (gx, gy), pts in grid.items():
        cand = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                cand.extend(grid.get((gx + dx, gy + dy), ()))
        for p in pts:
            best = float('inf')
            for q in cand:
                if q is p:
                    continue
                d = (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2
                if d < best:
                    best = d
            best = math.sqrt(best)
            if best < 0.02: o1 += 1
            if best < 0.1: o2 += 1
            if best < cell:
                nn_sum += best; nn_cnt += 1
    n = len(pos)
    print("%s %d %d %.1f %.1f %.4f %.4f %.3f %.4f %.4f %.4f" % (
        c.split('/')[-1], n, bad, max(xs) - min(xs), max(ys) - min(ys),
        pct(L, 0.5), pct(L, 0.9), sd / mean, o1 / n, o2 / n, nn_sum / max(nn_cnt, 1)))
