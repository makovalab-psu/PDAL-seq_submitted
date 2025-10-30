library(tidyverse)

state = c("0", "1", "2", "3", "4", "5", "6", "7","8","9")
nonb = c("all", "APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z")

df = expand.grid(state, nonb)
colnames(df) = c("State", "nonB")

Mean.observed = c()
Mean.null = c()
Wilcox.p = c()

for (i in 1:nrow(df)){
	File_path = paste("data/09_Human_HMM/04_Multivariate_guassian/10/State_nonB_enrichment/state_", df$State[i], "/", df$nonB[i], "_enrichment.bed", sep = "")
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

write.csv(df, "data/09_Human_HMM/04_Multivariate_guassian/10/State_nonB_enrichment/hsapiens.csv", row.names = F)
