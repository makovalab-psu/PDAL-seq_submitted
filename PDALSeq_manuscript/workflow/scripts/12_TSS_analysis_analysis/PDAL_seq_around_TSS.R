library(dplyr)
library(matrixStats)

args <- commandArgs(trailingOnly = TRUE)
#args = c("data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/1000_nucleotides/Hsapien_GM24385_PDALseq.bg", "data/12_TSS_analysis_analysis/scrambled/state_all/promoter_protCode_singleTrans.bed", "100")

df.bed = readr::read_delim(args[2], col_names = c("chr", "start", "end", "gene", "null", "strand"))

df.PDALseq = readr::read_delim(args[1], col_names = c("chr", "start", "end", "reads"))

Read_sum = sum(df.PDALseq$reads)/1000000
l.df.PDALseq = df.PDALseq %>% group_by(chr) %>% group_split()

chr_index = c()

for (i in 1:length(l.df.PDALseq)){
	chr_index[i] = l.df.PDALseq[[i]]$chr[1]
}

kb = as.numeric(args[3])
kb_range_plus = -kb:kb
kb_range_minus = kb:-kb

M = matrix(nrow = length(kb_range_plus), ncol = nrow(df.bed))

for (i in 1:nrow(df.bed)){
	#print(i)
	if (df.bed$strand[i] == "+"){ 
		kb_range = kb_range_plus
		start = df.bed$start[i]
	}
	if (df.bed$strand[i] == "-"){ 
		kb_range = kb_range_minus
		start = df.bed$end[i]
	}
	df.PDAL.seq.i = l.df.PDALseq[[which(chr_index == df.bed$chr[i])]]
	intervals = which.min(abs(df.PDAL.seq.i$start - start)) + kb_range
	intervals[which(intervals < 1)] = NA
	M[,i] = df.PDAL.seq.i$reads[intervals]
} 

Average_PDALseq = rowMeans(M, na.rm = T)/Read_sum
SD_PDALseq = rowSds(M, na.rm = T)/Read_sum
N = ncol(M)
kb_range = -kb:kb
df = data.frame(kb_range, Average_PDALseq, SD_PDALseq, N)

write.csv(df, args[4], row.names = F)
