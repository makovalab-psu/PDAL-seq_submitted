
for i in {0..9}; do
    total_length=$(awk '{sum += $3 - $2 + 1} END {print sum;}' data/09_Human_HMM/04_Multivariate_guassian/10/state_${i}.bed)
    echo "${i},$total_length"
done

