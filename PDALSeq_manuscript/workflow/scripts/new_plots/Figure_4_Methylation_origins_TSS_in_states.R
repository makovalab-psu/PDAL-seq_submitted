
library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)

####Define variables####

Exp.levels   = c("Hsapien_GM24385_PDALseq",
		 "Hsapien_GM24143_PDALseq",
		 "Hsapien_HEK293T_PDALseq",
		 "Hsapien_CCL86-Raji_PDALseq",
		 "Ptroglodytes_AG18354_PDALseq",
		 "Ptroglodytes_AG18358_PDALseq",
		 "Ppaniscus_KB8711_PDALseq",
		 "Ppaniscus_KB14571_PDALseq",
		 "Ggorilla_KB3781_PDALseq",
		 "Ggorilla_F7549_PDALseq",
		 "Ppygmaeus_AG05252_PDALseq",
		 "Pabelii_AG06213_PDALseq",
		 "Ssyndactylus_JAMBI_PDALseq",
		 "Ssyndactylus_Karenina_PDALseq",
		 "Hsapien_GM24385_TruSeq",
		 "Ptroglodytes_AG18354_TruSeq",
		 "Ppaniscus_KB8711_TruSeq",
		 "Ggorilla_KB3781_TruSeq",
		 "Ppygmaeus_AG05252_TruSeq",
		 "Pabelii_AG06213_TruSeq",
		 "Ssyndactylus_JAMBI_TruSeq")

Exp.label   = c("H. sapiens GM24385 PDAL-seq",
		 "H. sapiens GM24143 PDAL-seq",
		 "H. sapiens HEK293T PDAL-seq",
		 "H. sapiens CCL86-Raji PDAL-seq",
		 "P. troglodytes AG18354 PDAL-seq",
		 "P. troglodytes AG18358 PDAL-seq",
		 "P. paniscus KB8711 PDAL-seq",
		 "P. paniscus KB14571 PDAL-seq",
		 "G. gorilla KB3781 PDAL-seq",
		 "G. gorilla F7549 PDAL-seq",
		 "P. pygmaeus AG05252 PDAL-seq",
		 "P. abelii AG06213 PDAL-seq",
		 "S. syndactylus JAMBI PDAL-seq",
		 "S. syndactylus Karenina PDAL-seq",
		 "H. sapiens GM24385 TruSeq",
		 "P. troglodytes AG18354 TruSeq",
		 "P. paniscus KB8711 TruSeq",
		 "G. gorilla KB3781 TruSeq",
		 "P. pygmaeus AG05252 TruSeq",
		 "P. abelii AG06213 TruSeq",
		 "S. syndactylus JAMBI TruSeq")

NONB_levels = c("all", "APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS",  "GQ")
NONB_labels = c("all", "APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS",  "GQ")
n5mC_levels = c("0_0.2", "0.2_0.4", "0.4_0.6", "0.6_0.8", "0.8_1")
n5mC_labels = c("0 to 0.2", "0.2 to 0.4", "0.4 to 0.6", "0.6 to 0.8", "0.8 to 1")

####Hidden Markov model#####

#Plot that shows signal in each state

df = readr::read_csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/Average_PDALseq_in_GHMM_states.csv") %>% arrange(MEAN)
State.levels = unique(df$state)

df$state = factor(df$state, levels = State.levels)
df$Experiment = factor(df$Experiment, levels = Exp.levels, labels = Exp.label)

P.tile = ggplot(df, aes(x = Experiment, y = state, fill = MEAN)) +
	geom_tile() +
	scale_fill_gradient(low = "white", high = "red", name = "Mean PDAL-seq signal") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
	      axis.title.x = element_blank(),
	      legend.position = "top",
	      legend.title = element_text(size = 8),
	      legend.text = element_text(size = 6),
	      axis.title.y = element_blank(),
	      axis.text.y = element_blank()) 

####Make a plot that show methylation in each state####

df.5mC = read.csv("results/10_Methylation_analysis/Methylation_by_state_and_species.csv")
df.5mC$State = factor(df.5mC$State, levels = rev(State.levels))

df.5mC$Species = factor(df.5mC$Species, levels = rev(c("Hsapien","Ptroglodytes", "Ppaniscus",    "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus")), labels = rev(c("H. sapiens","P. troglodytes", "P. paniscus",    "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus")))

df.5mC$Range = factor(df.5mC$Range, levels = n5mC_levels, labels = n5mC_labels)

P.5mC = ggplot(df.5mC, aes(x = Range, y = Species, fill = Percent)) +
	facet_wrap(~State, ncol = 1, strip.position = "left") +
	geom_tile() +
	scale_fill_gradient(name =  "Percent (%)", low = "white", high = "black") +
	theme_classic() +
	scale_y_discrete(position = "right") +
	xlab("5mC rate") +
	theme(axis.text.x = element_text(angle = 90, size = 6),
	      axis.title.x = element_text(size = 8),
	      legend.position = "top",
	      legend.title = element_text(size = 8),
	      legend.text = element_text(size = 6),
	      strip.text = element_text(size = 6),
	      axis.title.y =  element_blank(),
	      axis.text.y = element_text(size = 6),
	      plot.title =  element_text(size = 8),
	      plot.margin = margin(t = 0, l = 0, r = 0, b = 0),
	      legend.key.size = unit(0.1, "in")) 
	      
df.hypo = df.5mC %>% filter(Range == "0 to 0.2")	


df.hypo$State = factor(df.hypo$State, levels = State.levels)

library(ggbeeswarm)

Species_palette = c("#9E9E9E", "#BE242F","#DC93B6", "#1C7D90", "#F2D983", "#FA9632")

P.hypo = ggplot(df.hypo , aes(x = Percent, y = State)) +
	geom_bar(stat="summary", fun = mean, color = "black", fill = "grey") +
    stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
    geom_beeswarm(mapping = aes(color = Species), size = 0.5) +
    theme_classic() +
    ggtitle("Hypo-methylated\nCpG sites") +
    xlab("Percent (%)") +
    scale_color_manual(values = rev(Species_palette)) +
    theme(axis.text.x = element_text(angle = 90, size = 6),
	      axis.title.x = element_text(size = 8),
	      legend.position = "none",
	      legend.title = element_text(size = 8),
	      legend.text = element_text(size = 6),
	      axis.title.y =  element_blank(),
	      axis.text.y =  element_blank(),
	      plot.title =  element_text(size = 8),
	      plot.margin = margin(t = 0.23, r = 0, b = 0.175, l = 0, unit = "in")) 
	      
P.A = plot_grid(P.5mC, P.hypo, rel_widths = c(3,1), nrow = 1)

 
####Look at signal around ORI####

df.OR = read.csv("results/11_Origin_analysis_analysis/PDALseq_signal_around_origins.csv")

df.OR$origin = factor(df.OR$origin, levels = c("core_origins_hglft_genome_e0596_1cc9e0", "stochastic_origins_hglft_genome_2ff9ed_1ccfd0"), labels = c("core origins\nn = 62892", "stochastic origins\nn = 254837"))
df.OR$state = factor(df.OR$state, levels = c("all", "4"), labels = c("all", "state 4"))

df.OR$sample = factor(df.OR$sample, levels = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq"), labels = c("GM24385-HG002", "GM24143-HG004", "HEK293T", "CCL86-Raji"))
df.OR$shuffle = factor(df.OR$shuffle, levels = c("unscrambled", "scrambled"))

df.OR = df.OR %>% filter(state != "state 4")

P.ORI = ggplot(df.OR, aes(x = kb_range, y = Average_PDALseq/1000, linetype = shuffle, color = origin)) +
	geom_line() +
	theme_classic() +
	facet_wrap(~sample , scales = "free") +
	scale_color_manual(values = viridis(3, option = "inferno")) +
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


P = plot_grid(P.A, P.ORI, nrow = 1, labels = c("A", "B"), label_size = 12)

ggsave("new_plots/Figure_4.svg", P, width = 6.5, height = 4, units = "in")

