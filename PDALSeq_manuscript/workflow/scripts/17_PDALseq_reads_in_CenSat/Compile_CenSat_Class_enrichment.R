library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

df.index = read.delim("PDALseq_samples.txt")

Censat = c("active_hor", "bSat", "cenSat", "ct", "dhor", "GAP", "GAP,HSat3", "gSat", "hor", "HSat1A", "HSat1B", "HSat2", "HSat2_3", "HSat3", "mixedSuperFamily", "mon", "mon-hor", "rDNA", "subTerm", "subTerm_mon", "subTerm_mon-hor")
Experiment = unique(df.index$Experiment)


df = expand.grid(Experiment, Censat)
colnames(df) = c("Experiment", "Censat")

l.Experiment = strsplit(as.character(df$Experiment), split = "_")

Biosample = c()
Species = c()
Technique = c()
Mean.observed = c()
Mean.null = c()
Wilcox.p = c()

for (i in 1:nrow(df)){
	File_path = paste("data/17_PDALseq_reads_in_CenSat/02_censat_enrichment_classes/", df$Censat[i], "/", df$Experiment[i], ".csv", sep = "")
	df.i = read.csv(File_path)
	if (length(unique(df.i$Observed)) == 1){
		Mean.observed[i] = mean(df.i$Observed)
		Mean.null[i] = mean(df.i$Null)
		Test = abs(Mean.observed[i] - Mean.null[i])
		P = length(which(abs(df.i$Null - Mean.null[i]) >= Test))/length(df.i$Null)
		Wilcox.p[i] = P
	}else{
		Mean.observed[i] = mean(df.i$Observed)
		Mean.null[i] = mean(df.i$Null)
		Wilcox = wilcox.test(df.i$Observed, df.i$Null)
		Wilcox.p[i] = Wilcox$p.value[[1]][1]
	}
	Biosample[i] = l.Experiment[[i]][2]
	Species[i] = l.Experiment[[i]][1]
	Technique[i] = l.Experiment[[i]][3]
}

df$Biosample = Biosample
df$Species = Species
df$Technique = Technique
df$Mean.observed = Mean.observed
df$Mean.null = Mean.null 
df$Wilcox.p = Wilcox.p
df$Enrichment = Mean.observed/Mean.null
df$log2Enrichment = log2(df$Enrichment)

write.csv(df, args[1], row.names = F)
