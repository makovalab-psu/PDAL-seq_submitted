#Filters a bam file into a new bam file depending on if the read maps inside a blacklist
#$1 The bam file you want to filter
#$2 The bam file you want the reads that map outside the black list to be placed
#$3 The bam file you want the reads that map outside the black list to be placed
#$4 The sample name used to find the blacklist
#$5 The taxon

Species=$(echo $4 | awk '{split($1, a, "_"); print a[1]}')

Blacklist=$(echo resources/alignments/$5/$Species.bed)

samtools view -b -L $Blacklist -U $2 $1 > $3

