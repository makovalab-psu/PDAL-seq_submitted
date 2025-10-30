line_number=$(grep -n -E 'chrX|chrY' $1 | head -n 1 | cut -d: -f1 | awk '{ print $1 - 1}')
echo $line_number
# Construct the sed command

head -n $line_number $2 > $3

