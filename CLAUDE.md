# LGL (Large Graph Layout) — thesis fork

Fork of Alex Adai's LGL, a force-directed layout tool for visualizing very
large graphs. This fork (copyright Fredrik Lindeberg, 2019/2020) was built to
generate Internet-topology graphs from BGP routing data for the PhD thesis
*"Coordinating the Internet: Thought Styles, Technology and Coordination"*
(Linköping University, 2021). Per the README, the repo's primary purpose is
to show the code used for the thesis for review — not to be a polished,
general-purpose library.

## Architecture

- **`src/`** — C++ core, built via autoconf/automake (`configure.ac`,
  `Makefile.am`). `lglayout.C` is the main layout engine: a multi-threaded
  (pthreads) particle-interaction physics simulation, depends on Boost.
  `lglbreakup.C`, `lglrebuild.C`, `lglfileconvert.C` handle graph
  pre/post-processing (splitting/rebuilding large graphs, format conversion).
- **`Java/`** — Java 2D viewer and image renderer. Produces `LGLLib.jar`,
  `LGLView.jar` (interactive 2D viewer), `ImageMaker.jar` (renders layout
  coordinates to large PNGs, up to 46340×46340px — limited by
  `java.awt.image.Raster`'s `Integer.MAX_VALUE` pixel-count cap). Built via
  `Java/build.sh`.
- **`perls/`** — Perl glue layer. `lgl.pl` is the driver script that calls
  the compiled C++ binaries; `LGLFormatHandler.pm` parses input edge files
  (ncol format: whitespace-delimited `<node1> <node2>` per line, plus other
  formats); `colorEdgesBasedOnLevel.pl` generates an edge-coloring file for
  the viewer based on hierarchy depth.
- **`scripts/`** — thesis-specific pipeline layered on top of generic LGL:
  downloads BGP RIB dumps (RIPE RIS `bview.*.gz` files), converts them with
  `parsebview.sh` into ncol graphs, then runs the full layout+render
  pipeline. Key entry points:
  - `creategraphfromdate.sh <year> <month>` — one-liner, fetches the last
    dump of that month from every RIS collector and produces a full
    graph+image.
  - `createseptemberseries.sh [first] [last]` — the same for September of
    every year (default 2000–now), deleting dumps once parsed (disk).
  - `creategraphfromurl.sh <url>` — same, from an arbitrary bview URL.
  - `create_run.sh <name>` → `bootstrap.sh` (parse+layout) → `run.sh`
    (drives `lgl.pl`, generates color file, renders PNG via `ImageMaker.jar`).
  - `parsebview.sh` parses each dump in parallel with `bin/bviewpairs`
    and sums the counts. `get_neigh.py`/`get_pairs.py` are the old
    bgpdump-based helpers it replaced; bviewpairs mirrors their semantics.
- **`tools/bviewpairs/`** — Rust (bgpkit-parser) tool: one MRT dump on
  stdin → AS–AS and prefix–AS pair counts. Same output as the old
  bgpdump pipeline except IPv6 prefixes are RFC 5952 formatted.
- **`setup.pl -i`** — orchestrates compiling the C++ and Java pieces and
  installing binaries into `./bin`.
- **`resources/images/`** — example rendered output (e.g.
  `internet_2016.png`) referenced from the README.
- **`testrun/`** — scratch output directories from past runs (per-run
  subfolders), not build artifacts to rely on.

## Build

```
make            # default target: Java jars only
make install    # jars, bviewpairs (cargo), C++ via setup.pl -i, all into ./bin
```

Requires: C++ compiler + Boost (Homebrew's is found automatically), Rust
(cargo), Perl 5+, Java 11+ (jars are built with `--release 11`), Python 3,
an X server only for the interactive viewers (ImageMaker runs headless).

CI (`.github/workflows/`) builds C++ (`make -C src all`), builds the Java
jars, then does a full `make install` + `creategraphfromdate.sh 2000 09` run
as an end-to-end smoke test, uploading the resulting PNG as an artifact.

## Notes for future work

- `run.sh` and `bootstrap.sh` have minimal error handling by design (per
  README) — don't assume a clean failure path.
- Coloring behavior lives in `perls/colorEdgesBasedOnLevel.pl`.
- Each graphing project needs its own directory under `testrun/` (a
  limitation inherited from the original LGL design), not a shared workspace.
- License: GPL (original LGL code, Alex Adai et al.); fork-specific changes
  also under GPL per LICENSE file.
