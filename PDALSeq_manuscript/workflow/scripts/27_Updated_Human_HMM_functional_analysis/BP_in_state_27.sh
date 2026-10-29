#!/bin/bash

for i in {0..9}; do
    total_length=$(awk '{sum += $3 - $2 + 1} END {print sum;}' \
        data/27_Updated_Human_HMM_functional_analysis/state_${i}.bed)
    echo "${i},$total_length"
done
