library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

Experiment = c("Hsapien_GM24385_PDALseq","Hsapien_GM24385_TruSeq","Ptroglodytes_AG18354_PDALseq","Ptroglodytes_AG18354_TruSeq","Ppaniscus_KB8711_PDALseq","Ppaniscus_KB8711_TruSeq","Ggorilla_KB3781_PDALseq","Ggorilla_KB3781_TruSeq","Pabelii_AG06213_PDALseq","Pabelii_AG06213_TruSeq","Ppygmaeus_AG05252_PDALseq","Ppygmaeus_AG05252_TruSeq","Ssyndactylus_JAMBI_PDALseq","Ssyndactylus_JAMBI_TruSeq")
df.censat = read.csv("data/17_PDALseq_reads_in_CenSat/00_scratch/Filtered_subclasses_for_analysis.csv")
Censat = unique(df.censat$file_path)


df = expand.grid(Experiment, Censat)
colnames(df) = c("Experiment", "Censat")

l.Experiment = strsplit(as.character(df$Experiment), split = "_")

Biosample = c()
Species = c()
Technique = c()
Class = c()
Subclass = c()
Mean.observed = c()
Mean.null = c()
Wilcox.p = c()

for (i in 1:nrow(df)){
	File_path = paste("data/17_PDALseq_reads_in_CenSat/04_censat_enrichment_subclasses/", df$Censat[i], "/", df$Experiment[i], ".csv", sep = "")
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
	v.Censat = strsplit(as.character(df$Censat[i]), split = "/")[[1]]
	Class[i] = v.Censat[1]
	Subclass[i] = v.Censat[2]
}

df$Biosample = Biosample
df$Species = Species
df$Technique = Technique
df$Class = Class 
df$Subclass = Subclass
df$Mean.observed = Mean.observed
df$Mean.null = Mean.null 
df$Wilcox.p = Wilcox.p
df$Enrichment = Mean.observed/Mean.null
df$log2Enrichment = log2(df$Enrichment)

write.table(df, args[1], sep = "\t", row.names = F, quote = F)
