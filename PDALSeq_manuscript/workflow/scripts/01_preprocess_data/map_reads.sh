#Maps reads to reference genome
#$1 the fastq file containing sequencing reads
#$2 the sample you are mapping, used to find the fasta
#$3 The output file
#$4 the number of threads to allocate
FASTA=$(echo $2 | awk '{split($1, a, "_"); print a[1]}')
echo resources/genomes/$FASTA.fasta

bwa mem -t $4 -o $3 resources/genomes/$FASTA.fasta $1
