state_bed=$1
annotation_bed=$2
output_file=$3
Subsamples=$4
Simulations=$5

temp_file="${output_file}_$$.temp"
shuffled_file="${output_file}_$$.shuffled"

echo "Observed,Null" > "${output_file}"

# Loop to run the process as many times as specified
for i in $(seq 1 "${Simulations}")
do
	# Randomly select lines
	shuf -n "${Subsamples}" "${annotation_bed}" > "${temp_file}"
	# Calculate Observed overlap
	Observed=$(bedtools intersect -a "${temp_file}" -b "${state_bed}" | awk '{sum += $4; n++} END {print sum;}')
	# Shuffle the BED file
	bedtools shuffle -i "${temp_file}" -g resources/genomes/Hsapien.chrom.sizes -chrom > "${shuffled_file}"
	# Calculate Null overlap
	Null=$(bedtools intersect -a "${shuffled_file}" -b "${state_bed}" | awk '{sum += $4; n++} END {print sum;}')
	# Save the results to the output file
	echo "$Observed,$Null" >> "${output_file}"
done

rm "${temp_file}"
rm "${shuffled_file}"

