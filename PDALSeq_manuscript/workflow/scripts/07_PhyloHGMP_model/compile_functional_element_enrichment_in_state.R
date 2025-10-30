library(tidyverse)

state = c("0", "1", "2", "3", "4", "5", "6", "7")
ann = c("cds_protCode_singleTrans", "ngnr_protCode_singleTrans", "chrG_cpgislands","promoter_protCode_singleTrans", "enhancer_protCode_singleTrans", "repeats_protCode_singleTrans", "exon_protCode_singleTrans", "rna_protCode_singleTrans", "gene_nonprotCode", "utr3_protCode_singleTrans", "gene_protCode_singleTrans", "utr5_protCode_singleTrans", "intron_protCode_singleTrans", "Recombination", "core_stochastic_origins_chm13v2.0_hglft.ucsc", "core_origins_hglft_genome_e0596_1cc9e0", "stochastic_origins_hglft_genome_2ff9ed_1ccfd0")

df = expand.grid(state, ann)
colnames(df) = c("State", "ann")

Mean.observed = c()
Mean.null = c()
Wilcox.p = c()

for (i in 1:nrow(df)){
	File_path = paste("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_functional_element_enrichment/state_", df$State[i], "/", df$ann[i], "_enrichment.bed", sep = "")
	df.i = read.csv(File_path)
	Mean.observed[i] = mean(df.i$Observed)
	Mean.null[i] = mean(df.i$Null)
	Wilcox = wilcox.test(df.i$Observed, df.i$Null)
	Wilcox.p[i] = Wilcox$p.value[[1]][1]
}

df$Mean.observed = Mean.observed
df$Mean.null = Mean.null 
df$Wilcox.p = Wilcox.p
df$Enrichment = Mean.observed/Mean.null
df$log2Enrichment = log2(df$Enrichment)

write.csv(df, "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_functional_element_enrichment/hsapiens.csv", row.names = F)

