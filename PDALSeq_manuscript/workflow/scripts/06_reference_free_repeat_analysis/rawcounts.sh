grep -v "@\|unit" $1 | awk '{if ($7>=75) {print;} }' | awk '{count[$3]+=1} END {for (word in count) print word, count[word]}' | sort > $2 #+1 for each occurence of an array
