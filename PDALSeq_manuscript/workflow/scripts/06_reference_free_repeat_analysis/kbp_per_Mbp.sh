#$1 #Original fasta file
#$2 #filtered data file containing repeats
#$3 # An output file with the columns repeat sequence, count, length, total Mbp in sam file, kbp repeat per Mbp sequenced, and original bam file  
bp=$(awk 'NR%2==0' $1 | wc -c)

echo $bp

cat $2 | awk -v bp="$bp" -v sample="$2" '{print $1 "\t" $2 "\t" bp "\t" sample}' > $3


