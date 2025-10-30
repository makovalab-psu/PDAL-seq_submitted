library(dplyr)
library(matrixStats)

args <- commandArgs(trailingOnly = TRUE)
#args = c("data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/1000_nucleotides/Hsapien_GM24385_PDALseq.bg", "resources/functional_annotations/Saswats_annotations/Hsapien/core_origins_hglft_genome_e0596_1cc9e0.bed", "100")

df.bed = readr::read_delim(args[2], col_names = c("chr", "start", "end", "label"))

df.PDALseq = readr::read_delim(args[1], col_names = c("chr", "start", "end", "reads"))
l.df.PDALseq = df.PDALseq %>% group_by(chr) %>% group_split()

chr_index = c()

for (i in 1:length(l.df.PDALseq)){
	chr_index[i] = l.df.PDALseq[[i]]$chr[1]
}

kb = as.numeric(args[3])
kb_range = -kb:kb

M = matrix(nrow = length(kb_range), ncol = nrow(df.bed))

for (i in 1:nrow(df.bed)){
	start = df.bed$start[i]
	df.PDAL.seq.i = l.df.PDALseq[[which(chr_index == df.bed$chr[i])]]
	intervals = which.min(abs(df.PDAL.seq.i$start - start)) + kb_range
	intervals[which(intervals < 1)] = NA
	M[,i] = df.PDAL.seq.i$reads[intervals]
} 

Average_PDALseq = rowMeans(M, na.rm = T) 
SD_PDALseq = rowSds(M, na.rm = T)
N = ncol(M)

df = data.frame(kb_range, Average_PDALseq, SD_PDALseq, N)

write.csv(df, args[4], row.names = F)
