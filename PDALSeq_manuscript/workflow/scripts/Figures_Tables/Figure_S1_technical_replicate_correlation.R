library(tidyverse)
library(viridis)

df = read.csv("results/03_Correlation_and_black_list_unmappable_regions/Technical_replicate_correlation.csv")



df.plot = df[-which(grepl("TruSeq", df$Sample_1)),]
df.plot = df.plot[-which(grepl("TruSeq", df.plot$Sample_2)),]

df.index = read.delim("PDALseq_samples.txt")

df.plot$Sample_1 = factor(df.plot$Sample_1, levels = df.index$Sequencing_sample)
df.plot$Sample_2 = factor(df.plot$Sample_2, levels = df.index$Sequencing_sample)


P = ggplot(df.plot, aes(x = Sample_1, y = Sample_2, fill = Spearman)) +
	facet_wrap(~Window) +	
	geom_tile() +
	theme_classic() +
	theme(axis.text.x = element_text(angle =90),
	      axis.text.y = element_blank(),
	      axis.title = element_blank()) +
	scale_fill_viridis(limits = c(0,1), option = "turbo") +
	coord_fixed()
	
ggsave("plots/Figure_S1_technical_replicate_correlation/Figure_S1.svg", P, width = 5, height = 10, scale = 3)

#fastq-dump --accession SRR11032812 --split-files --split-spot --gzip
