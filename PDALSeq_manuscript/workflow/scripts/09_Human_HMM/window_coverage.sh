#Calculates genome coverage from a genome windows file
#$1 The bam file with mapped reads
#$2 The windows file
#$3 The output file
#$4 The experiment

split -l 500 $2 temp/$4.ann.

Annotations=($(ls temp/$4.ann.*))

for i in "${Annotations[@]}"
do
   megadepth $1 --threads 1 --annotation ./$i --op sum > $i.reads
done

cat temp/$4.ann.*.reads | sort -k 1,1 -k2,2n > $3

rm temp/$4.ann.*
