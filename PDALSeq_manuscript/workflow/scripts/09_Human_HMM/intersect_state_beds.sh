#!/bin/bash
#$1 query bed file containing the human HMM state that you want to examine
#$2 database bed file contianing the phylogenetic HMM state you want to intersect 1 with

Total=$(awk '{sum+=$3-$2+1;} END{print sum;}' $1)

Intersect=$(bedtools intersect -a $1 -b $2 | awk '{sum+=$3-$2+1;} END{print sum;}')

echo $Total","$Intersect

