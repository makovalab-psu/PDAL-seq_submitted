#Calculates genome coverage from a genome windows file
#$1 The bam file with mapped reads
#$2 The taxon in the alignment
#$3 The output
#$4 The sample used to determine which annotation file to use

Species=$(echo $4 | awk '{split($1, a, "_"); print a[1]}')

Windows=$(echo resources/alignments/$2/$Species.bed)

mkdir temp/$2

split -l 500 $Windows temp/$2/$4.ann.

Annotations=($(ls temp/$2/$4.ann.*))

for i in "${Annotations[@]}"
do
   megadepth $1 --threads 1 --annotation ./$i --op sum > $i.reads
done

cat temp/$2/$4.ann.*.reads > $3

rm temp/$2/$4.ann.*
