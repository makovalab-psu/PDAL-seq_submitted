#Converts raw sequencing reads in the compressed fastq.gz format into raw sequencing reads in the fasta format
#$1 a fastq.gz input file path
#$2 a fasta output file path

gunzip -c $1 | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $2
