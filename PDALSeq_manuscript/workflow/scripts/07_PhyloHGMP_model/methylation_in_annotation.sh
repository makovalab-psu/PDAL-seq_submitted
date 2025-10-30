#$1 the annotattion file (bed)
#$2 The methylation file (bg)
#$3 The output file where the first column will be the number of Cs and the second will be the rate of methylation
bedtools intersect -a $2 -b $1 | awk '{ $4 = sprintf("%.3f", $4); print $4}' | sort | uniq -c > $3
