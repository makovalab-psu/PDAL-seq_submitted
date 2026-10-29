library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

K562_experiments = c("Hsapien_K562_PDALseq")

Exp = c(
	"Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq",
	"Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq",
	K562_experiments
)
Ann = c(
	"cds_protCode_singleTrans", "ngnr_protCode_singleTrans", "chrG_cpgislands",
	"promoter_protCode_singleTrans", "enhancer_protCode_singleTrans", "repeats_protCode_singleTrans",
	"exon_protCode_singleTrans", "rna_protCode_singleTrans", "gene_nonprotCode",
	"utr3_protCode_singleTrans", "gene_protCode_singleTrans", "utr5_protCode_singleTrans",
	"intron_protCode_singleTrans", "Recombination",
	"core_stochastic_origins_chm13v2.0_hglft.ucsc",
	"core_origins_hglft_genome_e0596_1cc9e0", "stochastic_origins_hglft_genome_2ff9ed_1ccfd0"
)

df = expand.grid(Exp, Ann)
colnames(df) = c("Exp", "Ann")

get_path = function(experiment, annotation){
	if (experiment %in% K562_experiments){
		paste("data/21_add_K562_to_analysis/12_read_enrichment/Read_enrichment/", experiment, "/", annotation, ".csv", sep = "")
	} else {
		paste("data/12_PDALseq/Read_enrichment/", experiment, "/", annotation, ".csv", sep = "")
	}
}

Mean.observed = c()
Mean.null = c()
Wilcox.p = c()

for (i in 1:nrow(df)){
	df.i = read.csv(get_path(df$Exp[i], df$Ann[i]))
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
