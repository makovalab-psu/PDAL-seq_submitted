
library(tidyverse)

read.size = function(x = "resources/genomes/Hsapien.chrom.sizes"){
	df = read.delim(x, header = F)
	return(sum(df$V2))
}

bed.size = function(x = "results/03_Correlation_and_black_list_unmappable_regions/blacklist_files/Hsapien.bed"){
	df = read.delim(x, header = F)
	return(sum(df$V3 - df$V2 + 1))
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
Blacklist.length = c()
Unique.maps.blacklist = c()
Unique.maps.outside.blacklist = c()

for (i in 1:nrow(df)){
	Total.length[i] = read.size(paste("resources/genomes/", df$Species[i], ".chrom.sizes", sep = ""))
	Blacklist.length[i] = bed.size(paste("results/03_Correlation_and_black_list_unmappable_regions/blacklist_files/", df$Species[i], ".bed", sep = ""))
	Unique.maps.blacklist[i] = calcUMR(paste("data/03_Correlation_and_black_list_unmappable_regions/08_estimate_unique_maps_by_blacklist/inside_blacklist/", df$Experiment[i], ".txt", sep = ""))
	Unique.maps.outside.blacklist[i] = calcUMR(paste("data/03_Correlation_and_black_list_unmappable_regions/08_estimate_unique_maps_by_blacklist/outside_blacklist/", df$Experiment[i], ".txt", sep = ""))
}

df$Total.length = Total.length
df$Blacklist.length = Blacklist.length
df$Unique.maps.blacklist = Unique.maps.blacklist
df$Unique.maps.outside.blacklist = Unique.maps.outside.blacklist


write.csv(df, "results/03_Correlation_and_black_list_unmappable_regions/Unique_maps_inside_and_outside_blacklist.csv", row.names = F)

