
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

State = c("state_0", "state_1", "state_2", "state_3", "state_4", "state_5", "state_6", "state_7", "state_8", "state_9")
Ann = c("Psuedogene", "LINE", "SINE", "LTR", "DNA_transposons", "Satelitte", "Simple_repeat", "Low_complexity")

df = expand.grid(State, Ann)
colnames(df) = c("State", "Ann")

Mean.observed = c()
Mean.null = c()
Wilcox.p = c()

for (i in 1:nrow(df)){
	File_path = paste("data/14_Repeat_masker_class_enrichment/02_annotations_enrichment/", df$State[i], "/", df$Ann[i], "_enrichment.bed", sep = "")
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
}

df$Mean.observed = Mean.observed
df$Mean.null = Mean.null 
df$Wilcox.p = Wilcox.p
df$Enrichment = Mean.observed/Mean.null
df$log2Enrichment = log2(df$Enrichment)

write.csv(df, args[1], row.names = F)
