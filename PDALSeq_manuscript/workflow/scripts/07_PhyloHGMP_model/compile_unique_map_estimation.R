
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

read.size = function(x = "resources/genomes/Hsapien.chrom.sizes"){
	df = read.delim(x, header = F)
	return(sum(df$V2))
}


bed.size = function(x = "resources/alignments/hominoidea/Hsapien.bed"){
	df = data.frame(readr::read_delim(x, col_names = c("chr", "start", "end", "strand", "chr.length", "sequence"), delim = "\t"))
	return(sum(df$end - df$start + 1))
}

calcUMR = function(x = "data/03_Correlation_and_black_list_unmappable_regions/08_estimate_unique_maps_by_blacklist/outside_blacklist/Hsapien_GM24385_PDALseq.txt"){
    df = read.delim(x, header = T,sep = "")
    HighMQ = which(df$Primary.score.minus.suplemental.score >= 10)
    if (length(HighMQ) > 0){
    	uniquely.mapped = sum(df$count[HighMQ])
    }else{
    	uniquely.mapped = sum(df$count)
    }
    return(uniquely.mapped/sum(df$count))
}


df.index = read.delim("PDALseq_samples.txt")

Experiment = unique(df.index$Experiment)

l.Experiment = strsplit(Experiment, split = "_")

Species = c()
Bio_sample = c()
Technique = c()

for (i in 1:length(Experiment)){
	Species[i] = l.Experiment[[i]][1]
	Bio_sample[i] = l.Experiment[[i]][2]
	Technique[i] = l.Experiment[[i]][3]
}

df = data.frame(Experiment, Species, Bio_sample, Technique)

Total.length = c()
Alignment.length = c()
Unique.maps.alignment = c()
Unique.maps.outside.alignment = c()

for (i in 1:nrow(df)){
	Total.length[i] = read.size(paste("resources/genomes/", df$Species[i], ".chrom.sizes", sep = ""))
	Alignment.length[i] = bed.size(paste("resources/alignments/", args[1], "/", df$Species[i], ".bed", sep = ""))
	Unique.maps.alignment[i] = calcUMR(paste("data/07_PhyloHGMP_model/04_estimate_unique_maps_in_alignment/", args[1], "/inside_alignment/", df$Experiment[i], ".txt", sep = ""))
	Unique.maps.outside.alignment[i] = calcUMR(paste("data/07_PhyloHGMP_model/04_estimate_unique_maps_in_alignment/", args[1], "/outside_alignment/", df$Experiment[i], ".txt", sep = ""))
}

df$Total.length = Total.length
df$Blacklist.length = Alignment.length
df$Unique.maps.alignment = Unique.maps.alignment
df$Unique.maps.outside.alignment = Unique.maps.outside.alignment


write.csv(df, args[2], row.names = F)

