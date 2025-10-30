library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

#args = c("data/07_PhyloHGMP_model/01_alignment_coverage/hominoidae/Hsapien_GM24385_PDALseq.bg", "data/07_PhyloHGMP_model/02_Normalize_alignment_coverage/hominoidae/Hsapien_GM24385_PDALseq.csv")

read.r = function(x = "data/07_PhyloHGMP_model/01_alignment_coverage/hominoidea/Hsapien_GM24385_PDALseq.bg"){
	df = data.frame(readr::read_delim(x, col_names = c("chr", "start", "end", "read.count"), delim = "\t"))
	return(df)
}

read.bed = function(x = "resources/alignments/hominoidea/Hsapien.bed"){
	df = data.frame(readr::read_delim(x, col_names = c("chr", "start", "end", "strand", "chr.length", "sequence"), delim = "\t")) %>% select(chr, start, end)
	return(df)
}

args[1]
args[2]
Species = gsub(".bg", "", strsplit(strsplit(args[1], split = "/")[[1]][5], split = "_")[[1]][1])
Species
Taxon = strsplit(args[1], split = "/")[[1]][4]
Taxon
	
df.r = read.r(args[1])
df.bed = read.bed(paste("resources/alignments/", Taxon, "/", Species, ".bed", sep = ""))
	
#Ensure that the bed file was not reordered during quantification
Same_order = FALSE
v.after = paste(df.r$chr, df.r$start, df.r$end)
v.before = paste(df.bed$chr, df.bed$start, df.bed$end)
v.match = which(v.after == v.before)
if (length(v.match) == length(v.after)){
	Same_order = TRUE
}
	
print("The bed file was not reordered during quantification")
print(Same_order)
	
#Normalize and compile data#
if (Same_order){
	#Tidy_formatted_data_frame
	df.r$Read_sum_billion = sum(df.r$read.count)/10^9
	df.r$alignment_length_bp = df.r$end - df.r$start + 1
	df.r$Normalized_read_count_per_bp_per_billion_reads = df.r$read.count/df.r$alignment_length_bp/df.r$Read_sum_billion
	df.out = df.r	
	Mean = mean(df.out$Normalized_read_count_per_bp_per_billion_reads)
	SD = mean(df.out$Normalized_read_count_per_bp_per_billion_reads)
	df.out$Normalized_read_count_per_bp_per_billion_reads = (df.out$Normalized_read_count_per_bp_per_billion_reads - Mean)/SD
}else{
	df.r$after = v.after
	l.df.r = df.r %>% group_by(chr) %>% group_split()
	chr.names = c()
	for (j in 1:length(l.df.r)){
		chr.names[j] = l.df.r[[j]]$chr[1]
	}
	read.count = rep(NA, nrow(df.bed))
	pb = txtProgressBar(min = 0, max = nrow(df.bed), initial = 0) 
	start_time <- Sys.time()
	for (j in 1:nrow(df.bed)){
		setTxtProgressBar(pb,j)
		df.j = l.df.r[[which(chr.names == df.bed$chr[j])]]
		read.count[j] = df.j$read.count[min(which(df.j$after == v.before[j]))]
	}
	end_time <- Sys.time()
	end_time - start_time
	df.bed$read.count = read.count
	df.bed$Read_sum_billion = sum(df.bed$read.count, na.rm = T)/10^9
	df.bed$alignment_length_bp = df.bed$end - df.bed$start + 1
	df.bed$Normalized_read_count_per_bp_per_billion_reads = df.bed$read.count/df.bed$alignment_length_bp/df.bed$Read_sum_billion
	Mean = mean(df.bed$Normalized_read_count_per_bp_per_billion_reads)
	SD = mean(df.bed$Normalized_read_count_per_bp_per_billion_reads)
	df.bed$Normalized_read_count_per_bp_per_billion_reads = (df.bed$Normalized_read_count_per_bp_per_billion_reads - Mean)/SD
	df.out = df.bed	
}	


write.csv(df.out, args[2], row.names = F)
