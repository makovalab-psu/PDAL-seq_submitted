"""Rename column 1 of a BED/bedGraph stream to UCSC native sequence names.

Usage: rename_chroms.py CHROM_ALIAS < in.bed > out.bed

CHROM_ALIAS is a UCSC chromAlias.txt file; its first column is the native name
and every other column is an alias. Lines on sequences with no match are dropped
and counted on stderr.
"""
import sys
from collections import Counter

alias = {}
with open(sys.argv[1]) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        names = line.rstrip("\n").split("\t")
        for name in names:
            if name:
                alias.setdefault(name, names[0])

dropped = Counter()
for line in sys.stdin:
    chrom, sep, rest = line.partition("\t")
    native = alias.get(chrom)
    if native is None:
        dropped[chrom] += 1
        continue
    sys.stdout.write(native + sep + rest)

for chrom, n in sorted(dropped.items()):
    sys.stderr.write(f"dropped {n} lines on {chrom} (not in chromAlias)\n")
