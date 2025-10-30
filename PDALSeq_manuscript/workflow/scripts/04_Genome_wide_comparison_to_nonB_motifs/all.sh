cat $1/* |
    sort -V -k1,1 -k2,2 |
    cut -f 1,2,3 |
    awk '{print $1 "\t" $2 "\t" $3}' > $2
