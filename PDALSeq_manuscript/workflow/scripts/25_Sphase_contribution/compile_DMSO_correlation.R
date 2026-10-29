library(tidyverse)

Samples = c("Hsapien_GM24385_LCL_PDALseq_js1014A20_1_S2_R1_001",
	"Hsapien_GM24385_LCL_PDALseq_js1014A20_2_S3_R1_001",
	"Hsapien_GM24385_LCL_PDALseq_js1014A40_1_S4_R1_001",
	"Hsapien_GM24385_LCL_PDALseq_js1014A40_2_S5_R1_001",
	"Hsapien_GM24143_LCL_PDALseq_js2007d20mM_S2_R1_001",
	"Hsapien_GM24143_LCL_PDALseq_js2007e20mM_S3_R1_001",
	"Hsapien_GM24143_LCL_PDALseq_js2007f40mM_S4_R1_001",
	"Hsapien_GM24143_LCL_PDALseq_js2007g40mM_S5_R1_001",
	"Hsa_GM24385_DMSO_PDALSeq_js4012a20mM_S1_R1_001",
	"Hsa_GM24385_DMSO_PDALSeq_js4012b20mM_S2_R1_001",
	"Hsa_GM24385_DMSO_PDALSeq_js4012c20mM_S3_R1_001",
	"Hsa_GM24143_DMSO_PDALSeq_js4012d20mM_S4_R1_001",
	"Hsa_GM24143_DMSO_PDALSeq_js4012e20mM_S5_R1_001",
	"Hsa_GM24143_DMSO_PDALSeq_js4012f20mM_S6_R1_001")

DMSO_samples = c("Hsa_GM24385_DMSO_PDALSeq_js4012a20mM_S1_R1_001",
	"Hsa_GM24385_DMSO_PDALSeq_js4012b20mM_S2_R1_001",
	"Hsa_GM24385_DMSO_PDALSeq_js4012c20mM_S3_R1_001",
	"Hsa_GM24143_DMSO_PDALSeq_js4012d20mM_S4_R1_001",
	"Hsa_GM24143_DMSO_PDALSeq_js4012e20mM_S5_R1_001",
	"Hsa_GM24143_DMSO_PDALSeq_js4012f20mM_S6_R1_001")

df = data.frame(expand_grid(Samples,Samples))
colnames(df) = c("Sample_1", "Sample_2")

Windows = c("1000000", "100000", "10000")

l_df = {}

for (i in 1:length(Windows)){
	l_df[[i]] = df
	l_df[[i]]$Window_size = Windows[i]
}

df = bind_rows(l_df)

read.bg = function(x){
	df = read.delim(x, sep = "", header = F)
	colnames(df) = c("chr", "start", "end", "reads")
	return(df)
}

get.path = function(sample, window){
	if (sample %in% DMSO_samples){
		paste0("data/25_Sphase_contribution/03_correlation/03_genome_window_coverage/",
		       window, "_nucleotides/", sample, ".bed")
	} else {
		paste0("data/03_Correlation_and_black_list_unmappable_regions/03_genome_window_coverage/",
		       window, "_nucleotides/", sample, ".bed")
	}
}

get.cell.line = function(sample){
	tokens = strsplit(sample, "_")[[1]]
	# DMSO samples: Hsa_GM24143_DMSO_... -> token index 2
	# Bulk samples: Hsapien_GM24143_LCL_... -> token index 2
	tokens[2]
}

Pearson  = rep(NA_real_, nrow(df))
Spearman = rep(NA_real_, nrow(df))

for (i in 1:nrow(df)){
	print(i)
	df1 = read.bg(get.path(df$Sample_1[i], df$Window_size[i]))
	df2 = read.bg(get.path(df$Sample_2[i], df$Window_size[i]))
	t.pearson  = cor.test(df1$reads, df2$reads, method = "pearson")
	Pearson[i] = t.pearson$estimate
	t.spearman  = cor.test(df1$reads, df2$reads, method = "spearman")
	Spearman[i] = t.spearman$estimate
}

df$Pearson  = Pearson
df$Spearman = Spearman

dir.create("results/25_Sphase_contribution", recursive = TRUE, showWarnings = FALSE)
write.csv(df, "results/25_Sphase_contribution/Technical_replicate_correlation.csv",
          row.names = FALSE)
   
