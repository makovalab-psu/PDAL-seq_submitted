
cat $1 | awk -v dir=$2 '{print $0 "\t" dir}' > $3
