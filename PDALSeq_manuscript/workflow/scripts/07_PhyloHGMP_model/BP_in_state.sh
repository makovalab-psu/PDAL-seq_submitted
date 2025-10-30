
for i in {0..7}; do
    total_length=$(awk '{sum += $3 - $2 + 1} END {print sum;}' data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_${i}.bed)
    echo "${i},$total_length"
done

