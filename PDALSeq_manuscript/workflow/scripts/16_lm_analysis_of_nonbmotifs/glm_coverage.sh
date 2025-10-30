#Calculates genome coverage from a genome windows file
#$1 The bam file with mapped reads
#$2 The output
#$3 The sample used to determine which annotation file to use

Species=$(echo $3 | awk '{split($1, a, "_"); print a[1]}')

Windows=$(echo data/16_lm_analysis_of_nonbmotifs/03_remove_sex_chromosomes_from_windows/1000_nucleotides/$Species.bed)

mkdir temp/$3_glm_coverage

split -l 10000 $Windows temp/$3_glm_coverage/Windows.ann.

Annotations=($(ls temp/$3_glm_coverage/Windows.ann.*))

for i in "${Annotations[@]}"
do
   megadepth $1 --threads 1 --annotation ./$i --op sum > $i.reads
done

cat temp/$3_glm_coverage/Windows.ann.*.reads > $2

rm -r temp/$3_glm_coverage
