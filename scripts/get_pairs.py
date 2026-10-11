#!/usr/bin/env python3

import sys

# python version of get_pairs

for line in sys.stdin:
    line = line.split()
    #line = line.split(' ')
    prefix = line[0]
    res = []
    for asn in line[1:]:
        if prefix == asn:
            continue
        else:
            res.append("{0} {1}".format(prefix, asn))

    if len(res) > 0:
        print("\n".join(res))

