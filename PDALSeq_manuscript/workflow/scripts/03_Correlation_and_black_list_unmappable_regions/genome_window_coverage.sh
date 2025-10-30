#Calculates genome coverage from a genome windows file
#$1 The bam file with mapped reads
#$2 The size of the windows file
#$3 The output
#$4 The sample

Species=$(echo $4 | awk '{split($1, a, "_"); print a[1]}')

Windows=$(echo data/03_Correlation_and_black_list_unmappable_regions/01_genome_windows/$2_nucleotides/$Species.bed)

echo $Windows

megadepth $1 --threads 1 --annotation ./$Windows --op sum > $3
