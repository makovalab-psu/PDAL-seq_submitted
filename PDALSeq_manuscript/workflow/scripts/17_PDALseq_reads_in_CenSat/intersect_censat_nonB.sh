#!/bin/bash
#$1 query bed file containing the censat annotation you want to examin
#$2 database bed file contianing the non-B motifs you want to calculate 

Total=$(bedtools merge -i $1 | awk '{sum+=$3-$2+1;} END{print sum;}')

Intersect=$(bedtools intersect -a $2 -b $1 | bedtools merge | awk '{sum+=$3-$2+1;} END{print sum;}')

echo $Total","$Intersect

