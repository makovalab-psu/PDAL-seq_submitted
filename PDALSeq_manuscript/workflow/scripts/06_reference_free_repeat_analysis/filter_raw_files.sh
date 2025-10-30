
file_counts=$1 #input file 1 is .rawcounts for filtering based presence of 100 repeats per 20 million reads
file_lengths=$2 #input file 2 is .rawlengths for filtering based on presence of 100 repeats per 20 million reads
output_file=$3 #output file 1 is column 1 = repeat, column 2 = count, and column 3 = length

#echo "counts file: " $file_counts
#echo "lengths file: " $file_lengths
#echo "output file: " $output_file

cat $file_lengths | awk '{print $2}' > $file_lengths.temp

#head $file_lengths.temp

pr -mts' ' $file_counts $file_lengths.temp > $output_file.temp

#head $output_file.temp

cat $output_file.temp | awk '{if ($2>=10) print;}' > $output_file

#head $output_file

rm $file_lengths.temp
rm $output_file.temp
