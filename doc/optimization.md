# Performance notes and remaining optimizations (2026-10)

Written during the LLM-assisted revival (Claude Code). Line numbers refer to the
tree as of this commit. Nothing in "Remaining" has been implemented.

## Where the time goes

A September series run (`scripts/createseptemberseries.sh`, 2000-2026, all RIS
collectors, full graph with prefixes) took ~7.5 h on a 12-core M-series Mac
before the layout optimizations below. Per year:

| Step | Share | Notes |
|---|---|---|
| `lglayout2D` (force simulation) | **94-96%** | 2025: 1.56M nodes, **44 levels**, 63 min at 10 threads with the old code, 12 min with the cell list |
| download + parse (`bviewpairs`) | 2-12 min | overlapped with the previous year's layout, so ~free |
| ImageMaker (3 PNGs at 2400x2400) | ~35 s | |
| `lglbreakup`, Perl conversion, `lglrebuild`, colour file | ~30 s | |

`lglayout2D` prints its own profile on stderr at the end of each pass: the time
per phase, the thread utilisation of the repulsion phase and the number of pair
tests. On the 2000 graph (108k nodes, 24 levels, 3800 iterations) with 10
threads, now 14.4 s wall (was 51 s before the cell list, 63 s at the start of the
revival):

| Phase | Time | Share |
|---|---|---|
| repulsion (cell-list pair kernel) | 9.4 s | 65%, thread utilisation 93% |
| pack (copy particles into the cell list) | 1.6 s | 11%, of which the serial `prepare` is a part |
| edges (attractive springs, atomic adds) | 1.75 s | 12% |
| integrate | 0.8 s | 5% |
| stats, placement, other | 0.25 s | 2% |

The pair kernel runs at 0.6-0.9 ns per pair test (NEON, four pairs at a time).
Visiting every pair from both sides costs 2.3x the pair tests of the old half
shell (12.9G vs 5.6G on a 25k-node test graph), which is still 3.5x faster than
the old scalar kernel at ~6 ns per test.

2025 (1.56M nodes, 44 levels, ~6000 iterations) with 10 threads: 705 s wall,
of which repulsion 561 s at 99.4% thread utilisation, pack 50 s (half of it the
sort by x), edges 43 s, integrate 33 s. It does 7.1 trillion pair tests (10.1
without the x window), i.e. ~750 per particle and iteration against ~175 on
the 2000 graph: the big hubs keep thousands of leaves in a few voxels for the
whole level, so 2025 is bound by the pair kernel alone and the small phases do
not matter there. On 2000 the x window is neutral (28% fewer lanes, but the
per-call setup weighs more with ~12 particles per voxel); on 2025 it gave -32%.

## Already done (don't redo)

- Mutex fix for macOS threads (`include/aPthread.hpp`), upstream Opte C++ changes
  (thread pool, atomic forces, `-L` getopt fix).
- Relaxed atomic `add2F`, `std::vector` voxel occupants, dynamic load balancing.
- **Cell list** (`include/cellList.hpp`): every iteration the particles of the
  voxels in the active box are copied into flat arrays (x, y, radius, anchor
  factor, pointer), one padded block per voxel. The repulsion
  (`calcInteractions`, `src/calcFuncs.C`) then takes work units of at most 256
  particles of one voxel, and for each particle sums the force from the 9 (2D)
  surrounding blocks with `repulsionFromBlock`: NEON on arm64 for 2D, a scalar
  loop otherwise. Every pair is seen from both sides, so a particle's force is
  summed by one thread and stored once (no atomics, no locks). Collisions are
  counted per lane and turned into noise from a per-thread generator
  (`ParticleInteractionHandler::noise`). Hubs with thousands of leaves in one
  voxel are split across threads by the units. Each block is sorted by x when
  packed, and for each particle only the window of a neighbouring block within
  the cutoff in x is visited (the window only moves forward, since the
  particles of the unit are sorted too). Blocks without anchors use a kernel
  without the radius, anchor-factor and ellipse loads when all radii are equal
  and there are no ellipse factors (`PairParams::plain`).
- Flat edge lists (`EdgeList` in `include/calcFuncs.h`, built per level by
  `buildEdgeLists`): the edge and stats phases take contiguous blocks instead of
  `std::advance` through boost's edge set per thread; the stats phase only sees
  the current level's edges. The integration phase takes contiguous vertex blocks
  and no lock (`clampIntegrateAndResetForce`).
- Parsing: `tools/bviewpairs` (Rust), one dump per core.
- Removed: `gridSchedual.hpp`, `voxelInteractionHandler.hpp` and the voxel lists
  per thread (the old static scheduling, which also allocated per-thread lists
  the size of the whole grid).

Behaviour differences from the thesis-era code that are not just summation
order: the collision noise comes from a per-thread xorshift generator instead of
the shared `std::rand()`; two colliding anchors no longer get noise (anchors never
drift); the cutoff test is done in float instead of double (the spring force is
zero at the cutoff distance, so a pair at the boundary contributes nothing
either way).

## Remaining, ranked by expected gain

Output quality must not get worse. "Equivalent" = only float summation order,
thread assignment or the random stream changes (layouts are stochastic anyway);
"changes output" needs a visual A/B against the thesis images in
`resources/images/latest/`.

1. **Iterations per level (changes output; biggest lever).** Every level runs until
   `iterationCtr > 150` (`beginSimulation`, `src/calcFuncs.C`); the convergence
   test compares a float sum of edge lengths against 1e-5/1e-6 and is mostly
   noise, so it rarely fires. 2025 has 44 levels; in 2013 the last 7 levels added
   2.4% of the nodes but cost ~40% of the work. The final "settle" pass
   (`src/lglayout.C`) is another ~6%. The per-level lines that `lglayout2D` now
   prints (nodes, iterations, seconds) show this directly. Options: scale the cap
   by the number of new nodes, merge small tail levels, sum dx in double with a
   looser tolerance.
2. **Fewer pair tests (equivalent).** The full shell still tests each pair twice
   (the x window removes the out-of-range part of the left and right columns
   only). A half shell would need per-thread force accumulators for the partner
   block or atomics again. Inside a dense hub every pair is within the cutoff, so
   there the n^2 is the physics. Measured and rejected: 1/sqrt via
   `vrsqrteq_f32` + 2 Newton steps instead of sqrt + div (no gain, the kernel is
   not division-bound), a mask-first pass (with 35% of lanes in range almost
   every group of 4 has one). Untested: two independent accumulator chains (8
   lanes per step), a y window as well by splitting blocks into y strips.
3. **Edge phase (equivalent, ~12%).** `onlyEdgeInteractions` does two atomic
   fetch_adds per edge (a CAS loop each on `boost::atomic<float>`). A per-vertex
   adjacency (CSR built with the edge list) would let each thread sum its own
   vertices' edge forces and store them once, like the repulsion does.
4. **Pack phase (equivalent, ~11%).** `CellList::prepare` walks the box serially
   (one occupancy read per box voxel, ~6k for 2000 but ~90k for 2025) and `fill`
   copies every particle. A parallel prefix sum, or keeping the SoA incrementally
   updated from `shift_particle`, would shrink it.
5. **Phase dispatch (equivalent).** Five thread-pool dispatches per iteration with
   a mutex and condition variable each; on the 2000 graph an iteration is ~4 ms.
   A persistent worker loop with a barrier would cut the idle time, which is now
   5-7% of the repulsion phase and more of the small phases.
6. **Reproducible layouts (equivalent).** With the cell list the forces no longer
   depend on which thread adds first, except through the voxel occupant order
   (`shift_particle` swaps) and the thread-bound noise generator. Sorting each
   block by particle index when packing and seeding the noise from the particle
   index and iteration would make a layout reproducible for any thread count,
   which the thesis-era code was not.
7. **Particle layout.** `Particle` is 160 B with an unused 72-byte `Amutex` and a
   duplicate `std::string id_`. Only the pack and integrate phases touch it now,
   so this moved down the list.
8. **Orchestration.** Two layouts in parallel with 5-6 threads each instead of one
   with 10; start the next layout while ImageMaker runs; skip `lgl.pl` part III
   (rebuild) when there is one connected set. With the layout 3.5x faster the
   download and parse (2-12 min per year) are no longer hidden behind it for the
   small years, so those matter more now.

Rejected: Barnes-Hut/FMM and multilevel layouts (FM3, sfdp) produce pictures that
are not comparable with the thesis; `-x` seeding from the previous year (10x faster
but clearly worse layouts, ~40% of nodes are new each early year); `-ffast-math`;
inlining the whole pair force into the voxel loops (slower than letting the
compiler decide, measured).

## Known bugs

Changing these changes the images compared with the thesis:
- `-L` never works as intended: `doesVertexHaveAnyChildren` (`src/calcFuncs.C`)
  looks only at the current level's edges, so it always returns false.
- `generatePlacementVector` (`src/calcFuncs.C`) scales "unit" vectors by their
  own length.

Do not affect current images:
- Downloads and parsing fail silently (no `wget` exit check or timeout, no
  `pipefail`); `run.sh` does not notice an `lgl.pl` failure (`$?` is tee's).
- `lglrebuild.C:136` reads `molecules[0]` after `pop_back` (UB); `lglbreakup.C`
  `(int) .5 * edgeCount` is always 0.
- ImageMaker drops edges beyond 10000 px, never draws vertices, and colours
  overlapping edges in a non-deterministic order.

## How to measure

- Build: `make install` (the 2D binary is `bin/lglayout2D`; `setup.pl` sets
  `DIMENSION` in `include/configs.h` per build). For one binary only:
  set `DIMENSION = 2` in `include/configs.h`, `make -C src clean lglayout`, and
  restore `configs.h`.
- Inputs: `testrun/bview_<year>_09/tmp/*/0.lgl` after a run (the series script
  deletes them, so copy one first). 2000 (108k nodes, ~15 s) for fast iteration,
  2025 (1.56M nodes) for the real case.
- Run: `bin/lglayout2D -o out.coords -y -O -t 10 -E 1x1 -e -l -L 0.lgl` (the flags
  `lgl.pl` uses). The profile is printed on stderr after the progress lines
  (`tr '\b' '\n' < err | grep -v ^Iteration`).
- Time: interleave old/new runs (A B A B ...), at least 3 each, report medians of
  wall and user+sys CPU. Keep the machine otherwise idle: on Apple Silicon a
  loaded machine moves single-thread runs to efficiency cores.
- Refactors that should not change the arithmetic: run old and new with `-t 1`
  on a small graph (a BFS subgraph of 25k nodes of 2000 takes ~35 s with the
  old code) and `cmp` the coords; single-threaded runs are deterministic.
- Quality: `python3 tools/layout_quality.py 0.lgl run1.coords run2.coords ...`
  prints nodes, NaN/Inf, extent, median/p90 edge length, CV of edge length and the
  share of nodes with a neighbour closer than 0.02/0.1. The new code must stay
  within the old code's run-to-run spread. Also render a PNG
  (`perl perls/colorEdgesBasedOnLevel.pl out.coords.edge_levels > colors`, then
  `java -Djava.awt.headless=true -jar Java/jar/ImageMaker.jar 2400 2400 0.lgl
  out.coords -c colors`, written next to the coords file) and compare by eye.
- Parsing changes must give byte-identical ncol files (`cmp`).
- Threads: build with `-fsanitize=thread` (also in `LIBS`) and run with
  `TSAN_OPTIONS=suppressions=file` where the file holds
  `race:boost::atomics::detail::non_atomic_load`; at `-O1` that false positive
  (the CAS loop of `boost::atomic<float>::fetch_add` starts with a plain load)
  is reported inside `onlyEdgeInteractions` instead and has to be read past.
