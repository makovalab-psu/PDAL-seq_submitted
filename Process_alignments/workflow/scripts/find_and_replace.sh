
cat $1 > $3

# Read the replacements file line by line
while IFS=',' read -r old new; do
    # Use sed to replace the old substring with the new substring in target.txt
    command="sed -e "\""s/$old/$new/g"\"" $3 > $3.tmp" 
    echo $command 
    eval $command 
    mv "$3.tmp" "$3"
done < $2

