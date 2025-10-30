####Look at TSS signal ####

df.TSS = read.csv("results/12_TSS_analysis_analysis/PDALseq_signal_around_TSS.csv")

df.TSS$state = factor(df.TSS$state, levels = c("all", "4", "7", "3", "1", "6", "0", "5","2"))

df.TSS$sample = factor(df.TSS$sample, levels = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq"), labels = c("GM24385", "GM24143", "HEK293T", "CCL86-Raji"))
df.TSS$shuffle = factor(df.TSS$shuffle, levels = c("unscrambled", "scrambled"))


P.TSS = ggplot(df.TSS, aes(x = kb_range - 1, y = Average_PDALseq, linetype = shuffle, color = sample)) +
	geom_line() +
	theme_classic() +
	facet_wrap(~state, scales = "free") +
	scale_color_manual(values = viridis(5)) +
	xlab("Distance from TSS (kb)") +
	ylab("PDAL-seq reads/(Million reads total)") +
	xlim(-25, 25)

ggsave("plots/Figure_3_Phylogenetic_PDALseq/TSS.svg", P.TSS, height = 6, width = 14)

