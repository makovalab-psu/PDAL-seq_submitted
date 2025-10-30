library(tidyverse)
state = c("0", "1", "2", "3", "4", "5", "6", "7")

df = expand.grid(state)
colnames(df) = c("State")

Mean.observed = c()
Mean.null = c()
Wilcox.p = c()

for (i in 1:nrow(df)){
	File_path = paste("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_RNA_expression/state_", df$State[i], ".txt", sep = "")
	df.i = read.csv(File_path)
	Mean.observed[i] = mean(df.i$Observed)
	Mean.null[i] = mean(df.i$Null)
	Wilcox = wilcox.test(df.i$Observed, df.i$Null)
	Wilcox.p[i] = Wilcox$p.value[[1]][1]
}

df$Mean.observed = Mean.observed
df$Mean.null = Mean.null 
df$Wilcox.p = Wilcox.p
df$Enrichment = Mean.observed/Mean.null
df$log2Enrichment = log2(df$Enrichment)


write.csv(df, "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_RNA_expression/hsapiens.csv", row.names = F)
