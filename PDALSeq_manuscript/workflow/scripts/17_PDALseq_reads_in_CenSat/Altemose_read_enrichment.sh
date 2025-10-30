state_bam=$1
output_file=$2
Subsamples=$3
Simulations=$4

Species=$(echo $5 | awk '{split($1, a, "_"); print a[1]}')
annotation_bed=$(echo data/17_PDALseq_reads_in_CenSat/07_Altemose_paper_censat_classes/$Species/$6.bed)
windows_bed=$(echo data/03_Correlation_and_black_list_unmappable_regions/01_genome_windows/1000_nucleotides/$Species.bed)

echo "Observed,Null" > $output_file

# Loop to run the process as many times as specified
for i in $(seq 1 $Simulations)
do
	# Check if annotation_bed is empty
    	if [ ! -s $annotation_bed ]; then
        	echo "," >> $output_file
        	continue  # Skip the rest of the loop iteration
    	fi
	
	# Randomly select lines
	bedtools intersect -a $windows_bed -b $annotation_bed | shuf -n $Subsamples | bedtools sort > $output_file.temp
	
	# Check if $output_file.temp is empty
    	if [ ! -s $output_file.temp ]; then
        	echo "," >> $output_file
        	continue  # Skip the rest of the loop iteration
    	fi
	
	# Calculate Observed overlap
	Observed=$(megadepth $state_bam --threads 1 --annotation ./$output_file.temp --op sum | awk '{sum += $4}END{print sum;}')
	# Shuffle the BED file
	bedtools shuffle -i $output_file.temp -g resources/genomes/$Species.chrom.sizes -chrom | bedtools sort > $output_file.shuffled
	# Calculate Null overlap
	Null=$(megadepth $state_bam --threads 1 --annotation ./$output_file.shuffled --op sum | awk '{sum += $4}END{print sum;}')
	# Save the results to the output file
	echo "$Observed,$Null" >> $output_file
	
	rm $output_file.temp
	rm $output_file.shuffled
done



