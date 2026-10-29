library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

K562_experiments = c("Hsapien_K562_PDALseq")

sample = c(
	"Hsapien_CCL86-Raji_PDALseq", "Hsapien_GM24385_PDALseq",
	"Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq",
	K562_experiments
)
shuffle = c("unscrambled", "scrambled")
state = c("all", "0", "1", "2", "3", "4", "5", "6", "7")
TSS = c("promoter_protCode_singleTrans")

df = expand.grid(sample, shuffle, state, TSS)
colnames(df) = c("sample", "shuffle", "state", "TSS")

get_path = function(sample, shuffle, state, tss){
	if (sample %in% K562_experiments){
		paste("data/21_add_K562_to_analysis/12_TSS/Mean_signal/", sample, "/", shuffle, "/state_", state, "/", tss, ".csv", sep = "")
	} else {
		paste("data/12_TSS_analysis_analysis/Mean_signal/", sample, "/", shuffle, "/state_", state, "/", tss, ".csv", sep = "")
	}
}

l.df = {}
for (i in 1:nrow(df)){
	df.i = read.csv(get_path(df$sample[i], df$shuffle[i], df$state[i], df$TSS[i]))
	df.i$sample = df$sample[i]
	df.i$shuffle = df$shuffle[i]
	df.i$state = df$state[i]
	df.i$TSS = df$TSS[i]
	l.df[[i]] = df.i
}

df.out = bind_rows(l.df)
write.csv(df.out, args[1], row.names = F)
