#All "State_4_all.bed"
bedtools merge -d 1 -i data/09_Human_HMM/04_Multivariate_guassian/10/state_7.bed > data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_all.bed

#Seperate out regions that intersect with ape states
bedtools intersect -a data/09_Human_HMM/04_Multivariate_guassian/10/state_7.bed -b data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_0.bed | bedtools sort | bedtools merge -d 1 > data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_0_intersect.bed
bedtools intersect -a data/09_Human_HMM/04_Multivariate_guassian/10/state_7.bed -b data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_1.bed | bedtools sort | bedtools merge -d 1 > data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_1_intersect.bed
bedtools intersect -a data/09_Human_HMM/04_Multivariate_guassian/10/state_7.bed -b data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_2.bed | bedtools sort | bedtools merge -d 1 > data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_2_intersect.bed
bedtools intersect -a data/09_Human_HMM/04_Multivariate_guassian/10/state_7.bed -b data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_3.bed | bedtools sort | bedtools merge -d 1 > data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_3_intersect.bed
bedtools intersect -a data/09_Human_HMM/04_Multivariate_guassian/10/state_7.bed -b data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_4.bed | bedtools sort | bedtools merge -d 1 > data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_4_intersect.bed
bedtools intersect -a data/09_Human_HMM/04_Multivariate_guassian/10/state_7.bed -b data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_5.bed | bedtools sort | bedtools merge -d 1 > data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_5_intersect.bed
bedtools intersect -a data/09_Human_HMM/04_Multivariate_guassian/10/state_7.bed -b data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_6.bed | bedtools sort | bedtools merge -d 1 > data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_6_intersect.bed
bedtools intersect -a data/09_Human_HMM/04_Multivariate_guassian/10/state_7.bed -b data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_7.bed | bedtools sort | bedtools merge -d 1 > data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_7_intersect.bed

#Regions where we expect high PDAL-seq

cat  data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_2_intersect.bed data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_3_intersect.bed data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_4_intersect.bed | bedtools sort | bedtools merge -d 1 > data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_high_PDAL-seq_expected.bed 

#Regions where PDAL-seq is variable between species

cat data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_1_intersect.bed data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_6_intersect.bed | bedtools sort | bedtools merge -d 1 > data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_variable_PDAL-seq_expected.bed 

#Regions where we do not expect high PDAL-seq base on other species

cat  data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_0_intersect.bed data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_7_intersect.bed data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_5_intersect.bed | bedtools sort | bedtools merge -d 1 > data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_low_PDAL-seq_expected.bed 

#Non syntenic reagions

cat data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_0_intersect.bed data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_1_intersect.bed data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_2_intersect.bed data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_3_intersect.bed data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_4_intersect.bed data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_5_intersect.bed data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_6_intersect.bed data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_Ape_state_7_intersect.bed | bedtools sort | bedtools merge -d 1 > data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_syntenic.bed 

bedtools intersect -v -a data/09_Human_HMM/04_Multivariate_guassian/10/state_7.bed -b data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_syntenic.bed | bedtools sort | bedtools merge -d 1 > data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/State_7_nonsyntenic.bed

