state_bed=$1
annotation_bed=$2
output_file=$3
Subsamples=$4
Simulations=$5

echo "Observed,Null" > $3

# Loop to run the process as many times as specified
for i in $(seq 1 $5)
do
	# Randomly select lines
	shuf -n $Subsamples $2 > $3.temp
	# Calculate Observed overlap
	Observed=$(bedtools intersect -a $3.temp -b $1 | awk '{sum += $4; n++} END {print sum;}')
	# Shuffle the BED file
	bedtools shuffle -i $3.temp -g resources/genomes/Hsapien.chrom.sizes -chrom > $3.shuffled
	# Calculate Null overlap
	Null=$(bedtools intersect -a $3.shuffled -b $1 | awk '{sum += $4; n++} END {print sum;}')
	# Save the results to the output file
	echo "$Observed,$Null" >> $3
done

rm $3.temp
rm $3.shuffled

