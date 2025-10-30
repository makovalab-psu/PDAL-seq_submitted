
library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)
library(tidyverse)
library(viridis)
library(cowplot)

Species = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Pabelii", "Ppygmaeus", "Ssyndactylus")
Species_labels = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. abelii", "P. pygmaeus", "S. syndactylus")
Species_palette = c("#9E9E9E", "#BE242F","#DC93B6", "#1C7D90", "#F2D983", "#FA9632", "#313757")


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



df = read.csv("results/17_PDALseq_reads_in_CenSat/SST1_SATR1_SATR2_LSAU_enrichment.csv")

df.VNTR_148 = read.delim("results/17_PDALseq_reads_in_CenSat/Censat_subclass_enrichment.txt") %>% filter(Subclass == "VNTR_148")

df = bind_rows(df, df.VNTR_148[,-c(6:7)])

df$Experiment = factor(df$Experiment, levels = Exp.levels, labels = Exp.labels)
df = df %>% filter(!is.na(Enrichment))
df$log2Enrichment[which(df$log2Enrichment > log2(3))] = log2(3)
df$log2Enrichment[which(df$log2Enrichment < -log2(3))] = -log2(3)
Sig = rep("ns", nrow(df))
Sig[which(df$Wilcox.p <= 0.05)] = "*"
Sig[which(df$Wilcox.p <= 0.01)] = "**"
df$Sig.labels = paste(round(df$Enrichment, digits = 1), Sig, sep = "\n")

df$Technique = factor(df$Technique, levels = c("PDALseq", "TruSeq"), labels = c("PDAL-Seq", "Control (TruSeq)"))

 P = ggplot(df, aes(x = Experiment, y = Censat, fill = log2Enrichment, label = Sig.labels )) +
	facet_wrap(~Technique, ncol = 1, scales = "free_x") +
	geom_tile() +
	geom_text(size = 2) +
	theme_classic() +
	scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-log2(3), log2(3))) + 
              theme(axis.text.x = element_text(size=6, color="black", angle=45, hjust = 1),
	      axis.title.x = element_blank(),
	      axis.text.y = element_text(size=6, color="black"),
	      axis.title.y = element_blank(),
	      plot.title = element_text(size=8, color="black"),
	      legend.key.size = unit(0.1, "in"),
    	      legend.position = "top",
	      legend.title = element_text(size=8, color="black"),
	      legend.text = element_text(size=6, color="black"),
	      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
	      
ggsave("new_plots/Figure_S18_select_satellite_subclasses.svg", P, width = 0.8*6.5, height = 0.8*8, units = "in")
	
