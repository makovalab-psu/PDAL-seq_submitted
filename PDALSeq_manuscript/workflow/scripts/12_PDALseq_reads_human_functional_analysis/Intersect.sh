state_bed=$1
annotation_bed=$2
output_file=$3


echo "Intersect_bp,Total_bp_state,Total_bp_ann" > $output_file

Total_state=$(cat $state_bed | awk '{sum += $3 - $2 + 1}END{print sum;}')
Total_ann=$(cat $annotation_bed | awk '{sum += $3 - $2 + 1}END{print sum;}')

# Calculate Observed overlap
Intersect=$(bedtools intersect -a $state_bed -b $annotation_bed | awk '{sum += $3 - $2 + 1}END{print sum;}')

# Save the results to the output file
echo "$Intersect,$Total_state,$Total_ann" >> $output_file


