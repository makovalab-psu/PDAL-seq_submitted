TOTAL_alignments=$(wc -l $1)
echo "Total_alignments= "$TOTAL_alignments
BP_Total=$(cat $1 | awk '{sum+=($3-$2+1)} END{print sum;}')
echo "BP_in_total_alignments= "$BP_Total

cat $1 | bedtools sort | bedtools merge > $1.tmp
Merged_alignments=$(wc -l $1.tmp)
echo "Merged_alignments= "$Merged_alignments
BP_merged=$(cat $1.tmp | awk '{sum+=($3-$2+1)} END{print sum;}')
echo "BP_in_merged_alignments= "$BP_merged

rm $1.tmp

BP_genome=$(cat $2 | awk '{sum+=($2)} END{print sum;}')
echo "BP_in_genome= "$BP_genome

