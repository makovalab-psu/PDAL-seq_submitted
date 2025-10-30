bedtools getfasta -tab -fi $2 -bed $1 | cut -f 2 | tr -d -c 'cCgG\n' | awk '{ print length; }' > $3
