library(tidyverse)

PDALseq = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq")

l.Normalized_PDALseq = {}

for (i in 1:length(PDALseq)){
	File_path = paste("data/09_Human_HMM/01_Remove_sex_chromosomes/", PDALseq[i], ".bg", sep = "")
	df.i = readr::read_delim(File_path, col_names = c("chr", "start", "end", "reads"))
	Read_sum_billion = sum(df.i$reads, na.rm = T)/10^9
	Length_bp = df.i$end - df.i$start + 1
	read_count_per_bp_per_billion_reads = df.i$reads/Length_bp/Read_sum_billion
	Mean = mean(read_count_per_bp_per_billion_reads)
	SD = mean(read_count_per_bp_per_billion_reads)
	Normalized_read_count = (read_count_per_bp_per_billion_reads - Mean)/SD
	l.Normalized_PDALseq[[i]] = Normalized_read_count
}

df.out = data.frame(l.Normalized_PDALseq)
colnames(df.out) = PDALseq

write.csv(df.out, "data/09_Human_HMM/02_Normalize_human_window_coverage/Normalized_human_window_coverage.csv", row.names = F)
