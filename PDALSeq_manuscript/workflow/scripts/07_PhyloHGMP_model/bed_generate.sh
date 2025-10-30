#$1 the repeat masker.out file
#$2 the class
#$3 The outfile

cat $1 | awk -v CLASS=$2 '$11 == CLASS {print $5 "\t" $6 "\t" $7}'
