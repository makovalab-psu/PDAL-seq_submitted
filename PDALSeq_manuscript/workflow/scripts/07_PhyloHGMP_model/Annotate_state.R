library(tidyverse)

options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
#args = c("resources/alignments/hominoidea/Hsapien.bed", "results/07_PhyloHGMP_model/alignment_clustering/hominoidea/10/GHMM_states.csv")

df.bed = readr::read_delim(args[1], delim = "\t", col_names = c("chr", "start", "end"))[,1:3]
df.state = readr::read_csv(args[2])

df.bed$hidden_state = df.state$hidden_state

l.df.bed = df.bed %>% group_by(hidden_state) %>% group_split()

for (i in 1:length(l.df.bed)){
	STATE = l.df.bed[[i]]$hidden_state[1]
	PATH = paste("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_", STATE, ".bed", sep = "")
	write.table(l.df.bed[[i]], PATH, sep = "\t", row.names = F, col.names = F, quote = F)
}
