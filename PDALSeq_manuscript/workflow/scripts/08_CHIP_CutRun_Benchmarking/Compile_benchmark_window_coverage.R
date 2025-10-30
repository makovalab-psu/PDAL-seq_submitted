library(dplyr)

options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)

#args = c("resources/Epigenetic_datasets/BG4_Rloop_CHIP_CutTag_experiments.txt", "results/08_CHIP_CutRun_Benchmarking/Genome_wide_coverage.csv")

df.index = read.csv(args[1], sep = "\t")

l.df = {}

for (i in 1:nrow(df.index)){
	file_path = paste("data/08_CHIP_CutRun_Benchmarking/06_genome_window_coverage/1000000_nucleotides/", df.index$Run[i], ".bg", sep = "")
	df.i = readr::read_delim(file_path,  col_names = c("chr", "start", "end", "reads"))
	Q99 = quantile(df.i$reads, 0.99)
	df.i$Normalized_reads = df.i$reads/Q99
	df.i$Normalized_reads[which(df.i$Normalized_reads > 1)] = 1
	df.i$Technique = df.index$Assay[i]
	df.i$Experiment = paste(df.index$Assay[i], df.index$Run[i])
	df.i$interval = 1:nrow(df.i)
	l.df[[i]] = df.i
}

df.HEK = readr::read_delim("data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/1000000_nucleotides/Hsapien_HEK293T_PDALseq.bg",  col_names = c("chr", "start", "end", "reads"))
Q99 = quantile(df.HEK$reads, 0.99)
df.HEK$Normalized_reads = df.HEK$reads/Q99
df.HEK$Normalized_reads[which(df.HEK$Normalized_reads > 1)] = 1
	
df.HEK$Technique = "PDAL-seq"
df.HEK$Experiment = "PDAL-seq"
df.HEK$interval = 1:nrow(df.HEK)

df = bind_rows(l.df, df.HEK)

write.csv(df, args[2], row.names = F)
