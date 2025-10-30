#Makes a bed file containing evenly spaced genome windows
#$1 a ucsc genome browser chrom sizes file
#$2 the output bed file
#$3 the window size in base pairs
bedtools makewindows -g $1 -w $3 | cut -f 1,2,3 | awk '{print $1 "\t" $2+1 "\t" $3}' > $2
