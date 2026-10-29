
library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)

####Look at TSS signal ####

df.TSS = read.csv("results/21_add_K562_to_analysis/PDALseq_signal_around_TSS.csv")

State.levels = rev(c("2", "5", "0", "6", "1", "3", "7", "4", "all")) 
State.labels = c("all", "8-origin", "7-recomb", "6-gene", "5-variable", "4-variable", "3-low", "2-low", "1-low")

df.TSS$state = factor(df.TSS$state, levels = State.levels, labels = State.labels)

df.TSS$sample = factor(df.TSS$sample, levels = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq", "Hsapien_K562_PDALseq"), labels = c("GM24385", "GM24143", "HEK293T", "CCL86-Raji", "K562"))
df.TSS$shuffle = factor(df.TSS$shuffle, levels = c("unscrambled", "scrambled"))


P.TSS = ggplot(df.TSS, aes(x = kb_range - 1, y = Average_PDALseq, linetype = shuffle, color = sample)) +
	geom_line() +
	theme_classic() +
	facet_wrap(~state) +
	scale_color_manual(values = viridis(6)) +
	xlab("Distance from TSS (kb)") +
	ylab("PDAL-seq reads/(Million reads total)") +
	xlim(-25, 25)

ggsave("revised_plots/Figure_S12_PDALSeq_around_human_TSS.svg", P.TSS, height = 4.5, width = 6.5)

