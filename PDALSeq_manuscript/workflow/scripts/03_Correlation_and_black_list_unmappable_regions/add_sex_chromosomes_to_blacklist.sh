cat $1 | grep -v "chrX" | grep -v "chrY" > $3.temp 
cat $2 | grep chrX | awk '{print $1 "\t0\t" $2}' > $3.chrX
cat $2 | grep chrY | awk '{print $1 "\t0\t" $2}' > $3.chrY

cat $3.temp $3.chrX $3.chrY > $3

rm $3.temp
rm $3.chrX
rm $3.chrY

