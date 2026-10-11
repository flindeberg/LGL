# Performance notes and remaining optimizations (2026-10)

Written during the LLM-assisted revival (Claude Code). Line numbers refer to the
tree as of this commit. Nothing in "Remaining" has been implemented.

## Where the time goes

A September series run (`scripts/createseptemberseries.sh`, 2000-2026, all RIS
collectors, full graph with prefixes) took ~7.5 h on a 12-core M-series Mac,
before the last two layout optimizations. Per year:

| Step | Share | Notes |
|---|---|---|
| `lglayout2D` (force simulation) | **94-96%** | 2025: 1.56M nodes, **44 levels**, ~46 min at 10 threads (was 63) |
| download + parse (`bviewpairs`) | 2-12 min | overlapped with the previous year's layout, so ~free |
| ImageMaker (3 PNGs at 2400x2400) | ~35 s | |
| `lglbreakup`, Perl conversion, `lglrebuild`, colour file | ~30 s | |

Inside `lglayout2D` (profile on the 2000 graph before the voxel-pair force sums
and the occupied-box walk): ~60% `calcInteractions` (voxel walk and pair checks),
~40% atomic force adds in the pair kernel, 3-4% edge loops and integration.
Re-profile before starting; the balance has shifted.

## Already done (don't redo)

- Mutex fix for macOS threads (`include/aPthread.hpp`), upstream Opte C++ changes
  (thread pool, atomic forces, `-L` getopt fix).
- Relaxed atomic `add2F` (`include/particle.hpp:112`), `std::vector` voxel occupants.
- Dynamic thread balancing: threads take chunks of voxels through an atomic counter.
- Forces summed per voxel pair before the atomic add
  (`include/voxelInteractionHandler.hpp:104`, `:129`;
  `include/particleInteractionHandler.hpp:207`).
- Only the voxels in the bounding box of placed particles are walked
  (`src/calcFuncs.C:34`, box built in `integrateParticles` at `:83`).
- Parsing: `tools/bviewpairs` (Rust), one dump per core.

## Remaining, ranked by expected gain

Output quality must not get worse. "Equivalent" = only float summation order or
thread assignment changes (layouts are stochastic anyway); "changes output" needs a
visual A/B against the thesis images in `resources/images/latest/`.

1. **Iterations per level (changes output; biggest lever).** Every level runs until
   `iterationCtr > 150` (`src/calcFuncs.C:290`); the convergence test compares a
   float sum of edge lengths against 1e-5/1e-6 and is mostly noise, so it rarely
   fires. 2025 has 44 levels; in 2013 the last 7 levels added 2.4% of the nodes but
   cost ~40% of the work. The final "settle" pass (`src/lglayout.C:348-350`) is
   another ~6%. Options: scale the cap by the number of new nodes, merge small tail
   levels, sum dx in double with a looser tolerance. Log iterations/dx per level
   first.
2. **Hub pile-ups (equivalent).** `-L` places all leaves of a parent within a tiny
   radius, so hubs with 4-6.5k leaves end up in one voxel: O(k²) pairs as one work
   item on one thread while the others wait. Split large voxels' pair loops into
   sub-tasks. Not measured; profile per-voxel time first.
3. **SoA cell list + branchless SIMD pair kernel (equivalent).** Occupants are
   pointers to ~152-byte `Particle`s; ~35% of distance checks pass, an
   unpredictable branch. Counting-sort particles by voxel into contiguous x/y
   arrays each iteration, run a masked NEON kernel, handle collisions/anchors in a
   scalar pass. Estimated 20-35%; high effort. Keep the double-precision `d² <
   cutoff` test so the set of interacting pairs is identical.
4. **Fixed-point int64 force atomics (equivalent, deterministic).** `f` is
   `boost::atomic<float>` (`include/particle.hpp:66`), whose `fetch_add` is a CAS
   loop. Scaled int64 with relaxed `fetch_add` compiles to `STADD` on arm64 and
   makes the sums order-independent. Smaller gain now that #voxel-pair sums exist.
5. **Edge phases (equivalent).** `onlyEdgeInteractions` (`src/calcFuncs.C:166`)
   and `collectEdgeStats` (`:127`) walk boost's edge list with `std::advance`, so
   every thread walks all E edges; `integrateParticles` does the same over
   vertices with interleaved (cache-line-sharing) assignment. Flatten to a
   `vector<pair<u32,u32>>` per level and give threads contiguous blocks. 2-5%.
6. **Particle layout (equivalent).** `Particle` is 152 B: a 72-byte `Amutex`
   (`include/particle.hpp:59`, only used by `enforceFLimit`, where each particle has
   one owner) and a duplicate `std::string id_` (`:60`). Hot fields (x, f, radius,
   anchor) span 2-3 cache lines. Shrink and reorder. 1-3%.
7. **Small things.** Per-thread RNG instead of `std::rand()` in `addNoise`
   (`include/particleInteractionHandler.hpp:151`, a data race on macOS and a global
   lock on glibc); compute each pair distance once instead of three times; PGO
   (5-15%, train on 2000).
8. **Orchestration.** Two layouts in parallel with 5-6 threads each instead of one
   with 10 (Amdahl estimate 1.2-1.4x throughput, but large layouts look
   memory-bound); start the next layout while ImageMaker runs; skip `lgl.pl` part III
   (rebuild) when there is one connected set.

Rejected: Barnes-Hut/FMM and multilevel layouts (FM3, sfdp) produce pictures that
are not comparable with the thesis; `-x` seeding from the previous year (10x faster
but clearly worse layouts, ~40% of nodes are new each early year); `-ffast-math`.

## Known bugs

Changing these changes the images compared with the thesis:
- `-L` never works as intended: `doesVertexHaveAnyChildren` (`src/calcFuncs.C:440`)
  looks only at the current level's edges, so it always returns false.
- `generatePlacementVector` (`src/calcFuncs.C:471`) scales "unit" vectors by their
  own length.
- Anchors can drift: noise is added to `p2` when both particles are anchors.

Do not affect current images:
- `gridSchedual.hpp:158` colouring is not conflict-free after the first pass (no
  longer used for scheduling).
- Downloads and parsing fail silently (no `wget` exit check or timeout, no
  `pipefail`); `run.sh` does not notice an `lgl.pl` failure (`$?` is tee's).
- `lglrebuild.C:136` reads `molecules[0]` after `pop_back` (UB); `lglbreakup.C`
  `(int) .5 * edgeCount` is always 0.
- ImageMaker drops edges beyond 10000 px, never draws vertices, and colours
  overlapping edges in a non-deterministic order.

## How to measure

- Build: `make install` (the 2D binary is `bin/lglayout2D`; `setup.pl` sets
  `DIMENSION` in `include/configs.h` per build).
- Inputs: `testrun/bview_<year>_09/tmp/*/0.lgl` after a series run; 2000 (108k
  nodes, ~45 s) for fast iteration, 2025 (1.56M nodes, ~45 min) for the real case.
- Run: `bin/lglayout2D -o out.coords -y -O -t 10 -E 1x1 -e -l -L 0.lgl` (the flags
  `lgl.pl` uses).
- Time: interleave old/new runs (A B A B ...), at least 3 each, report medians of
  wall and user+sys CPU. Keep the machine otherwise idle.
- Quality: `python3 tools/layout_quality.py 0.lgl run1.coords run2.coords ...`
  prints nodes, NaN/Inf, extent, median/p90 edge length, CV of edge length and the
  share of nodes with a neighbour closer than 0.02/0.1. The new code must stay
  within the old code's run-to-run spread. Also render a PNG
  (`perl perls/colorEdgesBasedOnLevel.pl out.coords.edge_levels > colors`, then
  `java -Djava.awt.headless=true -jar Java/jar/ImageMaker.jar 2400 2400 0.lgl
  out.coords -c colors`) and compare by eye.
- Parsing changes must give byte-identical ncol files (`cmp`).
