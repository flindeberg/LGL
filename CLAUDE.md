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
  `bgpdump` + `parsebview.sh` into ncol graphs, then runs the full
  layout+render pipeline. Key entry points:
  - `creategraphfromdate.sh <year> <month>` — one-liner, fetches a RIPE dump
    for that month and produces a full graph+image.
  - `creategraphfromurl.sh <url>` — same, from an arbitrary bview URL.
  - `create_run.sh <name>` → `bootstrap.sh` (parse+layout) → `run.sh`
    (drives `lgl.pl`, generates color file, renders PNG via `ImageMaker.jar`).
  - Python variants (`get_neigh.py`, `get_pairs.py`) are ~3x faster than
    their bash equivalents for the same parsing step.
- **`setup.pl -i`** — orchestrates compiling the C++ and Java pieces and
  installing binaries into `./bin`.
- **`resources/images/`** — example rendered output (e.g.
  `internet_2016.png`) referenced from the README.
- **`testrun/`** — scratch output directories from past runs (per-run
  subfolders), not build artifacts to rely on.

## Build

```
make            # builds Java jars + C++ binaries
make install    # also runs setup.pl -i, installs into $(PROJECTDIR)/bin
```

Requires: C++ compiler + Boost (version-sensitive — not fully future-proofed),
`bgpdump` (compile from source on Debian-likes, packaged version is often
broken), Perl 5+, Java 8+, Python 3, an X server for the graphical viewers.

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
