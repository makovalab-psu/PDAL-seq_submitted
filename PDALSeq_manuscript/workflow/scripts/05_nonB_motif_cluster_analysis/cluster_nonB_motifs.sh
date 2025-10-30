cat $1/* | sort -V -k1,1 -k2,2 | bedtools cluster -d $3 > $2

