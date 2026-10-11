#!/bin/bash


## Uses a bview file (such as http://data.ris.ripe.net/rrc00/latest-bview.gz),
## i.e. not updates but the entire view, and generates ncol files which can
## be used as graphs from that.

zcat=zcat # start with normal zcat on path

# if we have gzcat (gnu-version), lets use it
if hash gzcat 2>/dev/null; then
	zcat=gzcat
fi

## Get the current folder and use it for naming
folder=$(basename $(pwd))
asas=${folder}.as_as.ncol
preas=${folder}.prefix_as.ncol
export ncol=${folder}.full.ncol

## Quicker sort
export LC_ALL=C

## Each dump is parsed on its own, several in parallel (one per collector
## for a creategraphfromdate.sh run), and the per-dump counts are summed.
## Counts add up, so the result is the same as parsing all dumps as one
## stream, just faster.
##
## Per dump, bviewpairs (tools/bviewpairs, built by make) reads the MRT
## data and counts, in one pass:
##   AS -> AS: pairs of neighbouring ASes on each AS path
##   AS -> prefix: each prefix paired with its origin AS(es)
## This replaces bgpdump | cut | sed | get_neigh.py / get_pairs.py | awk,
## with the same output (IPv6 prefixes are written in RFC 5952 form).
##
## Then, for both: awk sums the per-dump counts and sort -n does a last
## numerical sort.

## make install puts bviewpairs in bin/, two levels up from the run directory
bviewpairs=${bviewpairs:-$(cd "$(dirname "$0")/../.." && pwd)/bin/bviewpairs}
if [ ! -x "$bviewpairs" ]; then
	echo "Error: $bviewpairs not found, run 'make install' first" >&2
	exit 1
fi

parse_one(){
	$zcat "$1" | "$bviewpairs" "$2.as_as" "$2.prefix_as"
}
export -f parse_one
export zcat bviewpairs

## One dump per core (gzcat and bviewpairs per dump, bviewpairs does the work)
jobs=${PARSE_JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu)}

parts=${folder}.parse_parts
rm -rf $parts && mkdir $parts

echo "in $(pwd)"
echo "Parsing $# dump(s) with $jobs in parallel for ASpath-info (i.e. AS -> AS routing) and prefix-info (i.e. AS -> prefix routing)"

i=0
time (for dump in "$@"; do
	i=$((i+1))
	printf '%s\0%s\0' "$dump" "$parts/$i"
done | xargs -0 -n 2 -P $jobs "$BASH" -c 'parse_one "$1" "$2"' _) ||\
	{ echo "Error: parsing failed, see above" >&2; rm -rf $parts; exit 1; }

echo "Summing per-dump counts into $asas and $preas"
time (cat $parts/*.as_as |\
	  awk '{c[$1 " " $2] += $3} END {for (k in c) print k " " c[k]}' | sort -n -S1G --parallel=24 > $asas
	cat $parts/*.prefix_as |\
	  awk '{c[$1 " " $2] += $3} END {for (k in c) print k " " c[k]}' | sort -n -S1G --parallel=24 > $preas)
rm -rf $parts

echo "Combining $asas and $preas into $ncol"
time (cat $asas $preas > $ncol)

