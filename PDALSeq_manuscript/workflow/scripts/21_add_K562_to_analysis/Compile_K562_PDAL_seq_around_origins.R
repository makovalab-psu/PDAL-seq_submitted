library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

K562_experiments = c("Hsapien_K562_PDALseq")

sample = c(
	"Hsapien_CCL86-Raji_PDALseq", "Hsapien_GM24385_PDALseq",
	"Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq",
	K562_experiments
)
shuffle = c("unscrambled", "scrambled")
state = c("all", "4")
origin = c("core_origins_hglft_genome_e0596_1cc9e0", "stochastic_origins_hglft_genome_2ff9ed_1ccfd0")

df = expand.grid(sample, shuffle, state, origin)
colnames(df) = c("sample", "shuffle", "state", "origin")

get_path = function(sample, shuffle, state, origin){
	if (sample %in% K562_experiments){
		paste("data/21_add_K562_to_analysis/11_origins/Mean_signal/", sample, "/", shuffle, "/state_", state, "/", origin, ".csv", sep = "")
	} else {
		paste("data/11_Origin_analysis_analysis/Mean_signal/", sample, "/", shuffle, "/state_", state, "/", origin, ".csv", sep = "")
	}
}

l.df = {}
for (i in 1:nrow(df)){
	df.i = read.csv(get_path(df$sample[i], df$shuffle[i], df$state[i], df$origin[i]))
	df.i$sample = df$sample[i]
	df.i$shuffle = df$shuffle[i]
	df.i$state = df$state[i]
	df.i$origin = df$origin[i]
	l.df[[i]] = df.i
}

df.out = bind_rows(l.df)
write.csv(df.out, args[1], row.names = F)
