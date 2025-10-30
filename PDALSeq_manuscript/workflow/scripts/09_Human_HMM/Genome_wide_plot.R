library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)


PDALseq = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq")

df.HMM = readr::read_csv("data/09_Human_HMM/03_Multivariate_guassian/8/GHMM_states.csv")

l.df = {}
l.P = {}

for (i in 1:length(PDALseq)){
	print(i)
	#bg_path = paste("data/09_Human_HMM/00_Normalize_human_window_coverage/", PDALseq[i], ".bg", sep = "")
	print("compiling data")
	bg_path = paste("data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/1000_nucleotides/", PDALseq[i], ".bg", sep = "")
	df.i = readr::read_delim(bg_path, col_names = c("chr", "start", "end", "reads"))
	df.i$interval = 1:nrow(df.i)
	df.i$HMM.reads = df.HMM[[i]]
	df.i$state = factor(as.character(df.HMM$hidden_state), levels = c("6", "0", "4", "1", "7", "5", "3", "2"))
	df.i$Experiment = PDALseq[i]
	Q99 = quantile(df.i$reads, 0.99)
	df.i$PDALseq = df.i$reads/Q99
	df.i$PDALseq[which(df.i$PDALseq > 1)] = 1
	print("Making plots")
	P.reads = ggplot(df.i, aes(x = interval, y = PDALseq, fill = state)) +
		ggtitle(PDALseq[i]) +
		geom_bar(stat = "identity") +
		theme_classic() +
		scale_fill_manual(values = viridis(8, option = "turbo")) 
	P.chr = ggplot(df.i, aes(x = interval, y = 1, fill = chr)) +
		geom_tile() +
		theme_void() +
		scale_fill_manual(values = viridis(length(unique(df.i$chr)))) +
		theme(legend.position = "none")
	
	print("Rendering plots")
	P = plot_grid(P.chr, P.reads, rel_heights = c(1,3), align = "v", ncol = 1)
	l.P[[i]] = P
	l.df[[i]] = df.i	
}

P.out = plot_grid(l.P[[1]], l.P[[2]], l.P[[3]], l.P[[4
]], ncol = 1, align = "v")

ggsave("test.png", P.out, height = 6, width = 10, scale = 2)

