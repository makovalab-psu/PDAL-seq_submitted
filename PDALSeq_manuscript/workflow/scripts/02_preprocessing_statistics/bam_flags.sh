#Counts the reads with different types of bam-flags
#$1 is a mapped bam file
#$2 is the output text file where the first column will be the number or reads and the second will be a sam flag
#Use https://broadinstitute.github.io/picard/explain-flags.html to interpret sam flags
samtools view $1 | cut -f 2 | sort | uniq -c > $2
