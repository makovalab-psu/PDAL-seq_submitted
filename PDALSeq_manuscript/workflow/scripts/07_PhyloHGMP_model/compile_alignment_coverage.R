
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

if (args[1] == "hominoidae"){
	TAXON = c("Hsapien","Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus")
}

if (args[1] == "hominidae"){
	TAXON = c("Hsapien","Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii")
}

if (args[1] == "homininae"){
	TAXON = c("Hsapien","Ptroglodytes", "Ppaniscus", "Ggorilla")
}

if (args[1] == "hominini"){
	TAXON = c("Hsapien","Ptroglodytes", "Ppaniscus")
}

df.index = read.delim("PDALseq_samples.txt") %>% filter(Species %in% TAXON)

Experiment = unique(df.index$Experiment)

l.Experiment = strsplit(Experiment, split = "_")

l.df = {}
l.M = {}

for (i in 1:length(l.Experiment)){
	df.i = readr::read_csv(paste("data/07_PhyloHGMP_model/02_Normalize_alignment_coverage/", args[1], "/", Experiment[i], ".csv", sep = ""))
	df.i$Species = l.Experiment[[i]][1]
	df.i$Bio_sample = l.Experiment[[i]][2]
	df.i$Technique = l.Experiment[[i]][3]
	df.i$Block = 1:nrow(df.i)
	df.i$Experiment = Experiment[i]
	l.df[[i]] = df.i
	M.i = matrix(df.i$Normalized_read_count_per_bp_per_billion_reads, nrow = 1)
	rownames(M.i) = Experiment[i]
	l.M[[i]] = M.i
}

df.out = bind_rows(l.df)

write.csv(df.out, args[2], row.names = F)

M.out = do.call(rbind, l.M) 

write.csv(t(M.out), args[3], row.names = F)
	
