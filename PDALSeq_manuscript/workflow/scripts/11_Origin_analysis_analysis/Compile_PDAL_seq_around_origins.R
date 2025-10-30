library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

sample = c("Hsapien_CCL86-Raji_PDALseq", "Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq")
shuffle = c("unscrambled", "scrambled")
state = c("all","4")
origin = c("core_origins_hglft_genome_e0596_1cc9e0", "stochastic_origins_hglft_genome_2ff9ed_1ccfd0")

df = expand.grid(sample, shuffle, state, origin)
colnames(df) = c("sample", "shuffle", "state", "origin")

l.df = {}
for (i in 1:nrow(df)){
	file_path = paste("data/11_Origin_analysis_analysis/Mean_signal/", df$sample[i], "/", df$shuffle[i], "/state_", df$state[i], "/", df$origin[i], ".csv", sep = "")
	df.i = read.csv(file_path)
	df.i$sample = df$sample[i]
	df.i$shuffle = df$shuffle[i]
	df.i$state = df$state[i]
	df.i$origin = df$origin[i]
	l.df[[i]] = df.i
}

df.out = bind_rows(l.df)
write.csv(df.out, args[1], row.names = F)

