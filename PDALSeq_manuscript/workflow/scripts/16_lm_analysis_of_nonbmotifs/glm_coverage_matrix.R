library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
#args = c("data/16_lm_analysis_of_nonbmotifs/05_glm_coverage/1000_nucleotides/Hsapien_GM24385_PDALseq.bg", "data/16_lm_analysis_of_nonbmotifs/07_glm_matrix/1000_nucleotides/Hsapien_GM24385_PDALseq.csv", "Hsapien_GM24385_PDALseq")

Species = strsplit(args[3], split = "_")[[1]][1]

####Read in the original bed####

bed_path = paste("data/16_lm_analysis_of_nonbmotifs/03_remove_sex_chromosomes_from_windows/1000_nucleotides/", Species, ".bed", sep = "")
df.bed = readr::read_delim(bed_path, delim = "\t", col_names = c("chr", "start", "end"))

Order_v = paste(df.bed$chr, df.bed$start, df.bed$end)

####Read in GC data####
gc_path = paste("data/16_lm_analysis_of_nonbmotifs/06_glm_GC_content/1000_nucleotides/", Species, ".txt", sep = "")
df.GC = readr::read_delim(gc_path, delim = "\t", col_names = c("GC.count"))

df.bed$GC_content = df.GC$GC.count/(df.bed$end - df.bed$start + 1)

####Read in the nonb data####

nonb = c("APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z")

for (i in 1:length(nonb)){
	#print(nonb[i])
	nonb_path = paste("data/16_lm_analysis_of_nonbmotifs/04_nonb_densisty_in_windows/1000_nucleotides/", Species, "/", nonb[i], ".bed", sep = "")
	df.nonb = readr::read_delim(nonb_path, delim = "\t", col_names = c("chr", "start", "end", "count", "length_nonb", "length_interval", "nonB_density"))
	df.bed[[(4+i)]] = df.nonb$nonB_density
	colnames(df.bed)[(4+i)] = nonb[i]
	#print(head(df.bed))
}

####Read in read count data####

read_path = args[1]

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






