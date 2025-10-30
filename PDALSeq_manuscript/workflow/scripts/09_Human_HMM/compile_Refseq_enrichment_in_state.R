library(tidyverse)

state = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
ann = c("meiotic_recombination_region", "mitotic_recombination_region", "non_allelic_homologous_recombination_region", "recombination_feature")

df = expand.grid(state, ann)
colnames(df) = c("State", "ann")

Mean.observed = c()
Mean.null = c()
Wilcox.p = c()

for (i in 1:nrow(df)){
	File_path = paste("data/09_Human_HMM/04_Multivariate_guassian/10/State_functional_refseq_enrichment/state_", df$State[i], "/", df$ann[i], "_enrichment.bed", sep = "")
	df.i = read.csv(File_path)
	print(head(df.i))
	Mean.observed[i] = df.i$Observed[1]
	Mean.null[i] = mean(df.i$Null)
	Xnull = abs(df.i$Null - Mean.null[i])
	Xobs = abs(Mean.observed[i] - Mean.null[i])
	Wilcox.p[i] = length(which(Xnull >= Xobs))/length(Xnull)
	print(Wilcox.p[i])
}

df$Mean.observed = Mean.observed
df$Mean.null = Mean.null 
df$Wilcox.p = Wilcox.p
df$Enrichment = Mean.observed/Mean.null
df$log2Enrichment = log2(df$Enrichment)

write.csv(df, "data/09_Human_HMM/04_Multivariate_guassian/10/State_functional_refseq_enrichment/hsapiens.csv", row.names = F)

