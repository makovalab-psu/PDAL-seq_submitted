
echo "count Primary.score.minus.suplemental.score"

samtools view $1 | \
    awk 'BEGIN {OFS="\t"}
        {
           for (i = 1; i <= NF; i++) {
                    if ($i ~ /^AS:i:/) as_value = substr($i, 6);
                    if ($i ~ /^XS:i:/) xs_value = substr($i, 6);
            }
          print as_value - xs_value
        }'  |
    sort |
    uniq -c
