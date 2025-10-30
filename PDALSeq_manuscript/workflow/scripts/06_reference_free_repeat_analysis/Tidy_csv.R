library(tidyverse)
library(viridis)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)

df.index = read.delim(args[1])

read.kbp_per_Mbp = function(x = 1){
	df.i.x = df.index[x,] 
	file_path = paste("data/06_reference_free_repeat_analysis/07_kbp_per_Mbp/", df.i.x$Sequencing_sample, ".dat", sep = "")
	df = read.delim(file_path, sep = "", header = F)
	if (ncol(df) == 6){colnames(df) = c("Repeat", "count", "test", "sequenced.bp", "test2", "File")}
	if (ncol(df) == 4){colnames(df) = c("Repeat", "count", "sequenced.bp", "File")}
	LENGTH = nchar(df$Repeat)
	df$kbp.per.mbp = 1000*df$count*LENGTH/df$sequenced.bp
	df$Sequencing_sample = df.i.x$Sequencing_sample
	df$Species = df.i.x$Species
	df$Bio_sample = df.i.x$Bio_sample
	df$Technique = df.i.x$Technique
	return(df)
}

l.df = lapply(1:nrow(df.index), read.kbp_per_Mbp)

df = bind_rows(l.df)

write.csv(df, args[2], row.names = F)


