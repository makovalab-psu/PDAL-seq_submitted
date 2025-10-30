library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

#args = c("Hsapien_GM24385_TruSeq", "data/05_nonB_motif_cluster_analysis/11_compile_nonb_cluster_coverage/200/Hsapien_GM24385_TruSeq.csv")

Species = strsplit(args[1], split = "_")[[1]][1]

####Read in the original bed####

bed_path = paste("data/05_nonB_motif_cluster_analysis/08_sex_chromosome_blacklist_removed_sorted_nonB_cluster_matrix/200/", Species, ".bed", sep = "")
df.bed = readr::read_delim(bed_path, delim = "\t", col_names = c("chr", "start", "end", "cluster.ID", "length.N", "MR", "TRI", "STR", "Z", "APR", "DR", "SLS", "GQ", "IR", "CRU"))
df.bed = df.bed[,1:3]

Order_v = paste(df.bed$chr, df.bed$start, df.bed$end)

####Read in GC data####
gc_path = paste("data/05_nonB_motif_cluster_analysis/10_nonB_cluster_matrix_GC_content/200/", Species, ".txt", sep = "")
df.GC = readr::read_delim(gc_path, delim = "\t", col_names = c("GC.count"))

df.bed$GC_content = df.GC$GC.count/(df.bed$end - df.bed$start + 1)

####Read in the nonb data####

nonb = c("APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z")

for (i in 1:length(nonb)){
	#print(nonb[i])
	nonb_path = paste("data/05_nonB_motif_cluster_analysis/11_nonb_densisty_in_clusters/200/", Species, "/", nonb[i], ".bg", sep = "")
	df.nonb = readr::read_delim(nonb_path, delim = "\t", col_names = c("chr", "start", "end", "cluster.ID", "length.N", "MR", "TRI", "STR", "Z", "APR", "DR", "SLS", "GQ", "IR", "CRU", "count", "length_nonb", "length_interval", "nonB_density"))
	df.bed[[(4+i)]] = df.nonb$nonB_density
	colnames(df.bed)[(4+i)] = nonb[i]
	#print(head(df.bed))
}

####Read in read count data####

read_path = paste("data/05_nonB_motif_cluster_analysis/09_nonb_cluster_coverage/200/", args[1], ".bg", sep = "")

df.reads = readr::read_delim(read_path, delim = "\t", col_names = c("chr", "start", "end", "read_count"))

Order_reads = paste(df.reads$chr, df.reads$start, df.reads$end)

Read_were_scrambled = T

if (length(which(Order_v == Order_reads)) == nrow(df.bed)){
	Read_were_scrambled = F
}

if (Read_were_scrambled){
	print("Reads were scrambled")
}else{
	df.bed$reads = df.reads$read_count
	write.csv(df.bed, file = args[2], row.names = F)
}



