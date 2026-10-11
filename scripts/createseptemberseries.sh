#!/bin/bash

## Creates graphs for the last BGP dump of September (i.e. the 30th) for
## every year in a range, using all RIS collectors and the full graph
## (AS + prefixes), the same way creategraphfromdate.sh does.
##
## Usage: ./createseptemberseries.sh [first_year] [last_year]
##        (defaults: 2000, the first September with RIS data, and this year)
##
## To keep disk usage down, each year's downloaded dumps are removed as soon
## as they have been parsed. The layout of one year runs in the background
## while the next year is downloaded and parsed.
##
## Images end up in each run directory and are also collected in
## ../testrun/september_series/

first=${1:-2000}
last=${2:-$(date +%Y)}

series=../testrun/september_series
mkdir -p $series

layout_pid=""

layout_year(){
    local y=$1
    local folder=bview_${y}_09
    cd ../testrun/$folder || return 1
    echo "[$(date '+%F %T')] $y: layout starting"
    ncol=${folder}.full.ncol ./run.sh -C > run_${folder}.log 2>&1
    if ! ls tmp/*/0.coords > /dev/null 2>&1; then
        echo "[$(date '+%F %T')] $y: layout FAILED, see testrun/$folder/run_${folder}.log"
        return 1
    fi
    # same as bootstrap.sh, minus the copy into resources/images
    cp tmp/*/0.coords*transparent.png ${folder}_transparent.png
    cp tmp/*/0.coords*light.png ${folder}_light.png
    cp tmp/*/0.coords*dark.png ${folder}_dark.png
    cp tmp/*/0.lgl ${folder}.lgl
    cp tmp/*/0.coords ${folder}.coords
    cp ${folder}_*.png ../september_series/
    ## drop intermediates that are only copies or inputs of lglbreakup,
    ## the large years would not fit on a small disk otherwise
    rm -f tmp/*_vertex_file_match tmp/*_new_lgl.lgl tmp/${folder}.full.lgl tmp/*/0.coords_2400x2400_*.png
    echo "[$(date '+%F %T')] $y: done ($(grep -a 'took' run_${folder}.log))"
}

for y in $(seq $first $last)
do
    folder=bview_${y}_09
    echo "[$(date '+%F %T')] $y: downloading and parsing"
    ## lgltest makes bootstrap.sh stop after parsing
    lgltest=1 ./creategraphfromdate.sh $y 09 > /dev/null 2>&1

    ## the dumps are no longer needed once parsed
    du -sh ../testrun/$folder/bview_${y}.09 2>/dev/null | sed "s/^/[$(date '+%F %T')] $y: removing dumps /"
    rm -f ../testrun/$folder/bview_${y}.09/bview_*.gz

    if [ ! -s ../testrun/$folder/${folder}.full.ncol ]; then
        echo "[$(date '+%F %T')] $y: parsing FAILED, see testrun/$folder/bootstrap_bview_${y}.09.log"
        continue
    fi
    echo "[$(date '+%F %T')] $y: parsed, $(wc -l < ../testrun/$folder/${folder}.full.ncol) edges"

    ## one layout at a time, it uses all cores
    if [ -n "$layout_pid" ]; then wait $layout_pid; fi
    (layout_year $y) &
    layout_pid=$!
done

if [ -n "$layout_pid" ]; then wait $layout_pid; fi
echo "[$(date '+%F %T')] all done, images in testrun/september_series/"
