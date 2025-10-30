library(tidyverse)
state = c("0", "1", "2", "3", "4", "5", "6", "7")

l.df = {}

for (i in 1:length(state)){
	file_path = paste("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/methylation/state_", state[i], ".txt", sep = "") 
	df.i = read.delim(file_path, header = F, sep = "")
	colnames(df.i) = c("count", "methylation_rate")
	df.i$percent = df.i$count/sum(df.i$count)
	df.i$state = state[i]
	l.df[[i]] = df.i
}

df = bind_rows(l.df)

write.csv(df, "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/methylation/hsapiens.csv", row.names = F)
