
chr=$(cat $1 | cut -f 1 | sort -T temp/ | uniq)

for x in $chr;
    do
        grep -w $x $1 > $2.$x.txt
	Rscript workflow/scripts/05_nonB_motif_cluster_analysis/nonB_cluster_matrix.R $2.$x.txt $2.$x.tmp
    done

cat $2.*.tmp > $2

rm $2.*.txt
rm $2.*.tmp 
