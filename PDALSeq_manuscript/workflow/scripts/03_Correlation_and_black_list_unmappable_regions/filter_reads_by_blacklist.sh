#Filters a bam file into a new bam file depending on if the read maps inside a blacklist
#$1 The bam file you want to filter
#$2 The bam file you want the reads that map outside the black list to be placed
#$3 The bam file you want the reads that map outside the black list to be placed
#$4 The sample name used to find the blacklist

Species=$(echo $4 | awk '{split($1, a, "_"); print a[1]}')

Blacklist=$(echo results/03_Correlation_and_black_list_unmappable_regions/blacklist_files/$Species.bed)

samtools view -b -L $Blacklist -U $2 $1 > $3

