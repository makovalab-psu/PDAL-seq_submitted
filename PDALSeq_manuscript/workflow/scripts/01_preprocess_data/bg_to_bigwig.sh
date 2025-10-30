#Maps reads to reference genome
#$1 the bedgraph file containing sequencing reads
#$2 the sample you are mapping, used to find the fasta
#$3 The output file

FASTA=$(echo $2 | awk '{split($1, a, "_"); print a[2]}' | awk '{split($1, a, "/"); print a[2]}')

echo $FASTA

bedGraphToBigWig $1 resources/genomes/$FASTA.chrom.sizes $3
