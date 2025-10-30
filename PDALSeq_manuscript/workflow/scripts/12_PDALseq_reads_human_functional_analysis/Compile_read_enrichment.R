library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

Exp = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq")
Ann = c("cds_protCode_singleTrans", "ngnr_protCode_singleTrans", "chrG_cpgislands","promoter_protCode_singleTrans", "enhancer_protCode_singleTrans", "repeats_protCode_singleTrans", "exon_protCode_singleTrans", "rna_protCode_singleTrans", "gene_nonprotCode", "utr3_protCode_singleTrans", "gene_protCode_singleTrans", "utr5_protCode_singleTrans", "intron_protCode_singleTrans", "Recombination", "core_stochastic_origins_chm13v2.0_hglft.ucsc", "core_origins_hglft_genome_e0596_1cc9e0", "stochastic_origins_hglft_genome_2ff9ed_1ccfd0")

df = expand.grid(Exp, Ann)
colnames(df) = c("Exp", "Ann")

Mean.observed = c()
Mean.null = c()
Wilcox.p = c()

for (i in 1:nrow(df)){
	File_path = paste("data/12_PDALseq/Read_enrichment/", df$Exp[i], "/", df$Ann[i], ".csv", sep = "")
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

write.csv(df, args[1], row.names = F)

