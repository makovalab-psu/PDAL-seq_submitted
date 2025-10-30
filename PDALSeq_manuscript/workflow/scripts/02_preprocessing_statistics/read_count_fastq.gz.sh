#Counts reads in a fastq.gz
#$1 the fastq.gz file that you want to count
#$2 the text file to store the anwer

echo $(zcat $1|wc -l)/4|bc > $2
