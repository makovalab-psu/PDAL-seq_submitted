#$1 input gff file
#$2 input feature to subset by
#$3 output bed file

grep $2 $1 | awk -F "\t" '{print $1 "\t" $4 "\t" $5}' > $3
