library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)

####Look at signal around ORI####

df.OR = readr::read_csv("results/22_Alternative_origins/PDALseq_signal_around_alternative_origins.csv")
unique(df.OR$sample)
unique(df.OR$origin)

origin.levels = c("GSM5658908_Ini-seq2.called.replication.origins",
                  "GSM6127777_Ini-seq2.called.replication.origin.15.min",
                  "GSM5658909_SNS.seq.fractions.0.5.2.kb_called.replication.origins",
                  "GSM5658910_SNS.seq.fractions.2.4.kb_called.replication.origins")
origin.labels = c("Ini-Seq2 GSM5658908",
                  "Ini-Seq2 GSM6127777",
                  "SNS-Seq GSM5658909",
                  "SNS-Seq GSM5658910")
                         
df.OR$origin = factor(df.OR$origin, levels = origin.levels, labels = origin.labels)
df.OR$state = factor(df.OR$state, levels = c("all", "4"), labels = c("all", "state 4"))

df.OR$sample = factor(df.OR$sample, levels = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq", "Hsapien_K562_PDALseq"), labels = c("GM24385-HG002", "GM24143-HG004", "HEK293T", "CCL86-Raji", "K562"))
df.OR$shuffle = factor(df.OR$shuffle, levels = c("unscrambled", "scrambled"))

df.OR = df.OR %>% filter(state == "all")

P.ORI = ggplot(df.OR, aes(x = kb_range, y = Average_PDALseq/1000, linetype = shuffle, color = origin)) +
	geom_line() +
	theme_classic() +
	facet_wrap(~sample , scales = "free", nrow = 1) +
	scale_color_manual(values = viridis(4, option = "inferno")) +
	xlab("Distance from ORI (kb)") +
	ylab("PDAL-Seq reads (103)") +
	guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) + 
    	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
	      axis.title.x = element_text(size = 8),
	      legend.position = "bottom",
	      legend.title = element_text(size = 8),
	      legend.text =  element_text(size = 6),
	      axis.title.y =  element_text(size = 8),
	      axis.text.y =  element_text(size = 6),
	      plot.title =  element_text(size = 8)) 
	      
ggsave("revised_plots/Figure_S11_PDALSeq_signal_around_alternative_origin_sets.svg", P.ORI, width = 6.5, height = 3.5, units = "in")

