library(tidyverse)

Windows = c("1000000", "100000", "10000")

K562_samples_for_corr = c(
	"Hsa_K562_LCL_PDALSeq_js4003_20mM1_S2_R1_001",
	"Hsa_K562_LCL_PDALSeq_js4003_20mM2_S3_R1_001",
	"Hsa_K562_LCL_PDALSeq_js4003_40mM1_S4_R1_001",
	"Hsa_K562_LCL_PDALSeq_js4003_40mM2_S5_R1_001"
)

read.bg = function(x){
	df = read.delim(x, sep = "", header = F)
	colnames(df) = c("chr", "start", "end", "reads")
	return(df)
}

get.path = function(sample, window){
	if (sample %in% K562_samples_for_corr){
		paste("data/21_add_K562_to_analysis/03_correlation/03_genome_window_coverage/", window, "_nucleotides/", sample, ".bed", sep = "")
	} else {
		paste("data/03_Correlation_and_black_list_unmappable_regions/03_genome_window_coverage/", window, "_nucleotides/", sample, ".bed", sep = "")
	}
}

df.index = read.delim("PDALseq_samples.txt")
Hsapien_existing = df.index %>% filter(Species == "Hsapien") %>% pull(Sequencing_sample) %>% as.character()

Sample_2_pool = c(K562_samples_for_corr, Hsapien_existing)

df.pairs = expand.grid(Sample_1 = K562_samples_for_corr, Sample_2 = Sample_2_pool, stringsAsFactors = FALSE)
df.pairs$Species = "Hsapien"
df.pairs$Comparison = ifelse(df.pairs$Sample_2 %in% K562_samples_for_corr, "K562_vs_K562", "K562_vs_Hsapien_other")

l.df.window = list()
for (i in 1:length(Windows)){
	df.i = df.pairs
	df.i$Window = Windows[i]
	l.df.window[[i]] = df.i
}

df = bind_rows(l.df.window)

Pearson = rep(NA, nrow(df))
Spearman = rep(NA, nrow(df))

for (i in 1:nrow(df)){
	print(i)
	Path1 = get.path(df$Sample_1[i], df$Window[i])
	Path2 = get.path(df$Sample_2[i], df$Window[i])
	df1 = read.bg(Path1)
	df2 = read.bg(Path2)
	t.pearson = cor.test(df1$reads, df2$reads, method = "pearson")
	Pearson[i] = t.pearson$estimate
	t.spearman = cor.test(df1$reads, df2$reads, method = "spearman")
	Spearman[i] = t.spearman$estimate
}

df$Pearson = Pearson
df$Spearman = Spearman

dir.create("results/21_add_K562_to_analysis", recursive = TRUE, showWarnings = FALSE)
write.csv(df, "results/21_add_K562_to_analysis/Technical_replicate_correlation.csv", row.names = F)
