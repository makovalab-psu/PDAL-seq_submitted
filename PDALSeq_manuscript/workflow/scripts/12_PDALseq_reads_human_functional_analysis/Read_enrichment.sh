state_bam=$1
annotation_bed=$2
output_file=$3
Subsamples=$4
Simulations=$5

echo "Observed,Null" > $output_file

# Loop to run the process as many times as specified
for i in $(seq 1 $Simulations)
do
	# Randomly select lines
	shuf -n $Subsamples $annotation_bed | bedtools sort > $output_file.temp
	# Calculate Observed overlap
	Observed=$(megadepth $state_bam --threads 1 --annotation ./$output_file.temp --op sum | awk '{sum += $4}END{print sum;}')
	# Shuffle the BED file
	bedtools shuffle -i $output_file.temp -g resources/genomes/Hsapien.chrom.sizes -chrom | bedtools sort > $output_file.shuffled
	# Calculate Null overlap
	Null=$(megadepth $state_bam --threads 1 --annotation ./$output_file.shuffled --op sum | awk '{sum += $4}END{print sum;}')
	# Save the results to the output file
	echo "$Observed,$Null" >> $output_file
done

rm $3.temp
rm $3.shuffled

