library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)
library(ggbeeswarm)

CenSat_levels = c("active_hor", "hor", "dhor", "mon-hor", "mon", "cenSat", "ct", "bSat", "gSat", "HSat1A", "HSat1B", "HSat2", "HSat2_3", "HSat3", "subTerm", "subTerm_mon", "subTerm_mon-hor", "mixedSuperFamily", "GAP", "GAP,HSat3", "GAP-HSat3", "rDNA")

nonb_levels = c("all", "APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS", "GQ")

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

Exp.labels   = c("H. sapiens GM24385 PDAL-seq",
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
		 
Species_levels = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Pabelii", "Ppygmaeus", "Ssyndactylus")	
Species_labels = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. abelii", "P. pygmaeus", "S. syndactylus")		 

df.enr = read.csv("results/17_PDALseq_reads_in_CenSat/Censat_class_enrichment.csv")
df.enr$Experiment = factor(df.enr$Experiment, levels = rev(Exp.levels), labels = rev(Exp.labels))
df.enr$Censat = factor(df.enr$Censat, levels = CenSat_levels )

df.int = read.csv("results/17_PDALseq_reads_in_CenSat/Censat_nonB_intersect.csv")
df.int$species = factor(df.int$species, levels = Species_levels, labels = Species_labels)
df.int$censat = factor(df.int$censat, levels = CenSat_levels )
df.int$Percent = 100*df.int$intersect/df.int$total
df.int$nonb = factor(df.int$nonb, levels = nonb_levels)

Exclude = c("GAP", "GAP,HSat3", "GAP-HSat3", "rDNA")
NonB_plot = c("APR", "IR", "TRI", "Z", "DR", "GQ")

####Summary/rDNA plot####

df.enr$rDNA = F
df.enr$rDNA[which(df.enr$Censat == "rDNA")] = T

P.compare = ggplot(df.enr, aes(x = Experiment, y = Enrichment)) +
	geom_bar(data = df.enr %>% filter(!rDNA), stat = "summary", fun = "mean") +
	stat_summary(data = df.enr %>% filter(!rDNA), geom = "errorbar", fun.data = mean_se, position = "dodge") +
	geom_jitter(mapping = aes(color = rDNA)) +
	scale_color_manual(values = c("black","red")) +
	theme_classic() +
	geom_hline(yintercept = 1, color = "blue") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
	      legend.position = "top",
	      axis.title.x = element_blank()) +
	scale_x_discrete(limits =Exp.labels) +
	scale_y_continuous(trans = "log2")

####Enrichment####


df.guide = df.enr %>% filter(!Censat %in% Exclude) %>% arrange(Censat)
x_lim = as.character(unique(df.guide$Censat))

Sig = rep("ns", nrow(df.enr))
Sig[which(df.enr$Wilcox.p <= 0.05)] = "*"
Sig[which(df.enr$Wilcox.p <= 0.01)] = "**"

Enr_labels = paste(round(df.enr$Enrichment, digits = 2), Sig)
Enr_labels[which(is.na(df.enr$Enrichment))] = "ns"
df.enr$Enr_labels = Enr_labels

P.enr = ggplot(df.enr %>% filter(!Censat %in% Exclude) %>% filter(Technique != "TruSeq"), aes(x = Censat, y = Experiment, fill = Enrichment, label = Enr_labels)) +
	geom_tile() +
	geom_text(color = "white") +
	scale_fill_viridis(option = "turbo") +
	theme_classic() +
	theme(axis.text.x = element_blank(),
	      legend.position = "right",
	      axis.title.x = element_blank()) +
	      scale_x_discrete(limits = x_lim)

####Intersect####

P.int = ggplot(df.int %>% filter(!censat %in% Exclude) %>% filter(nonb %in% NonB_plot), aes(x = censat, y = Percent, fill = species)) +
	facet_wrap(~nonb, ncol = 1, scales = "free_y", strip.position = "left") +
	geom_bar(stat = "identity", position = "dodge") +
	scale_fill_manual(values = viridis(7, option = "turbo")) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
	      panel.grid.major = element_line(color = "grey80"),
	      legend.position = "right",
	      axis.title.x = element_blank())+
	      scale_x_discrete(limits = x_lim)


P.grid = plot_grid(P.enr, P.int, ncol = 1, align = "v")

ggsave("plots/Figure_X_Read_enrichment_in_censat_classes/Overall_distribution_highlight rDNA.svg", P.compare, width = 6.5, height = 6, units = "in")
ggsave("plots/Figure_X_Read_enrichment_in_censat_classes/Enrichment_and_nonb_motifs.svg", P.grid, width = 13, height = 10, units = "in")
