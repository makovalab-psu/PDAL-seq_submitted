#Calculates genome coverage from a genome windows file
#$1 The bam file with mapped reads
#$2 The size of the windows file
#$3 The output
#$4 The sample used to determine which annotation file to use

Species=$(echo $4 | awk '{split($1, a, "_"); print a[1]}')

Windows=$(echo data/05_nonB_motif_cluster_analysis/08_sex_chromosome_blacklist_removed_sorted_nonB_cluster_matrix/$2/$Species.bed)

split -l 10000 $Windows temp/$4.ann.

Annotations=($(ls temp/$4.ann.*))

for i in "${Annotations[@]}"
do
   megadepth $1 --threads 1 --annotation ./$i --op sum > $i.reads
done

cat temp/$4.ann.*.reads | bedtools sort -faidx resources/genomes/$Species.chrom.sizes > $3

rm temp/$4.*
