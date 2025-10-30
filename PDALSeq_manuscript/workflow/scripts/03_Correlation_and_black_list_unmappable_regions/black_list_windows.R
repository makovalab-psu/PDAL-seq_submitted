library(tidyverse)
library(viridis)
library(cowplot)

args = commandArgs(trailingOnly = TRUE)
Threshold = args[1]

read.bed = function(x){
    df = read.delim(x, header = F)
    colnames(df) = c("chr", "start", "end", "count")
    return(df)
}

blacklist = function(x = "data/03_Correlation_and_black_list_unmappable_regions/04_genome_window_coverage_for_blacklist/High_MapQ/Hsapien_GM24385_TruSeq.bg", y = "data/03_Correlation_and_black_list_unmappable_regions/04_genome_window_coverage_for_blacklist/All_reads/Hsapien_GM24385_TruSeq.bg"){
	df.High_MapQ = read.bed(x)
	df.all = read.bed(y)

	df = df.all[,1:3]
	df$High_MapQ = df.High_MapQ$count
	df$All = df.all$count
	df$Ratio = df$All/df$High_MapQ
	df$Ratio[which(is.infinite(df$Ratio))] = 2
	df$Ratio[which(is.na(df$Ratio))] = 2
	df$Ratio[which(df$Ratio > 2)] = 2
	df = df[which(df$Ratio >= Threshold),]
	return(df)
}

Species = c( "Hsapien", "Ptroglodytes", "Ppaniscus",  "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus")

TruSeq = c("Hsapien_GM24385_TruSeq", "Ptroglodytes_AG18354_TruSeq", "Ppaniscus_KB8711_TruSeq", "Ggorilla_KB3781_TruSeq",  "Ppygmaeus_AG05252_TruSeq",    "Pabelii_AG06213_TruSeq",  "Ssyndactylus_JAMBI_TruSeq")

for (i in 1:length(Species)){
	print(i)
	high = paste("data/03_Correlation_and_black_list_unmappable_regions/04_genome_window_coverage_for_blacklist/High_MapQ/", TruSeq[i], ".bg", sep = "")
	all  = paste("data/03_Correlation_and_black_list_unmappable_regions/04_genome_window_coverage_for_blacklist/All_reads/", TruSeq[i], ".bg", sep = "")
	out  = paste("data/03_Correlation_and_black_list_unmappable_regions/05_raw_blacklist_files/", Species[i], ".bed", sep = "")
	df.bl = blacklist(high, all)
	write.table(df.bl, out, sep = "\t", row.names = F, col.names = F, quote = F)
}
