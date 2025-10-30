library(tidyverse)
library(viridis)
library(cowplot)

args = commandArgs(trailingOnly = TRUE)
#args = "plots/Figure_S5_CHIP_CutRun_Benchmarking/Figure_S5.svg"

df = readr::read_csv("results/08_CHIP_CutRun_Benchmarking/Genome_wide_coverage.csv")

Experiments = unique(df$Experiment)

df.index = data.frame(expand_grid(Experiments, Experiments))
colnames(df.index) = c("Sample_1", "Sample_2")

Spearman = rep(NA, nrow(df.index))

for (i in 1:nrow(df.index)){
	df.x = df %>% filter(Experiment == df.index$Sample_1[i])
	df.y = df %>% filter(Experiment == df.index$Sample_2[i])
	test = cor.test(df.x$Normalized_reads, df.y$Normalized_reads, method = "spearman")
	Spearman[i] = test$estimate
	print(paste(df.index$Sample_1[i], df.index$Sample_2[i], Spearman[i]))
}

df.index$Spearman = Spearman

Samples = unique(df.index$Sample_1)
Samples = Samples[rev(order(Samples))]
Samples = c(Samples[-which(Samples == "PDAL-seq")], "PDAL-seq")
df.index$Sample_1 = factor(df.index$Sample_1, levels = Samples)
df.index$Sample_2 = factor(df.index$Sample_2, levels = Samples)

P.cor = ggplot(df.index, aes(x = Sample_1, y = Sample_2, fill = Spearman, label = round(Spearman, digits = 2))) +
	geom_tile() +
	geom_text(color = "white") +
	scale_fill_viridis(option = "turbo") +
	theme_classic() +
	theme(axis.title = element_blank(),
	      axis.text.x = element_text(angle = 90, hjust = 1)) 

CHR = unique(df$chr)
df$chr = factor(df$chr, levels = CHR)
P = ggplot(df %>% filter(Experiment %in% c("PDAL-seq", "R-loop CUT&Tag SRR26244778", "G4 CUT&Tag SRR25010686")), aes(x = interval, y = Normalized_reads, fill = chr)) +
	geom_bar(stat = "identity") +
	theme_classic() +
	facet_wrap(~Experiment, ncol = 1) +
	scale_x_continuous(limits = c(min(df$interval), max(df$interval)), expand = c(0, 0)) +
	theme(legend.position = "bottom",
	      axis.text.x = element_blank()) +
	ylab("Normalized reads per Mbp")

P.out = plot_grid(P.cor, P, ncol = 1, labels = c("A", "B"))

ggsave(args[1], P.out, height = 7, width = 6, scale = 2, limitsize = FALSE)

