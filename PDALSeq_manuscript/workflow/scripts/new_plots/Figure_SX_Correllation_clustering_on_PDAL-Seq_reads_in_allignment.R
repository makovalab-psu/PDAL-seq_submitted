
library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)

####A-B Correllation plots####

df.cor = read.csv("data/07_PhyloHGMP_model/alignment_clustering/hominoidae/PDALseq_correllation_in_alignment.csv")

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
		 
df.cor$Sample_1 = factor(df.cor$Sample_1, levels = Exp.levels, labels = Exp.label)     
df.cor$Sample_2 = factor(df.cor$Sample_2, levels = Exp.levels, labels = Exp.label)

df.cor$Comparison = "PDAL-seq & TruSeq"
df.cor$Comparison[which(grepl("PDAL-seq", df.cor$Sample_1) & grepl("PDAL-seq", df.cor$Sample_2))] = "PDAL-seq & PDAL-seq"
df.cor$Comparison[which(grepl("TruSeq", df.cor$Sample_1) & grepl("TruSeq", df.cor$Sample_2)) ] = "TruSeq & TruSeq"

P.cor.top = ggplot(df.cor, aes(x = Sample_1, y = spearman, color = Comparison, shape = Comparison)) +
	geom_beeswarm() +
	theme_classic() +
	scale_y_continuous(limits = c(-0.5,1), breaks = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), labels = c("-0.5", "", "0", "", "0.5", "", "1")) +
	scale_color_manual(values = c("red", "blue", "grey")) +
	theme(axis.text.x = element_text(angle = 90, size =6),
    	  axis.text.y = element_text(size = 6),
    	  axis.title.x = element_blank(),
    	  axis.title.y = element_text(size = 8),
    	  strip.text = element_text(size = 8),
          legend.position = "bottom",
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8))



P.cor.bottom = ggplot(df.cor, aes(x = Sample_1, y = Sample_2, fill = spearman)) +
	geom_tile() +
	scale_fill_viridis(option = "turbo") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, size =6),
    	  axis.text.y = element_blank(),
    	  axis.title.x = element_text(size = 8),
    	  axis.title.y = element_text(size = 8),
    	  strip.text = element_text(size = 8),
          legend.position = "right",
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8)) +
          xlab("Sample 1") +
          ylab("Sample 2")
          
	     

####B PCA plot####

df.PCA = read.csv("data/07_PhyloHGMP_model/alignment_clustering/hominoidae/PDALseq_PCA_in_alignment.csv")

df.PCA$Species = factor(df.PCA$Species, levels = unique(df.PCA$Species))

library(ggrepel)

P.PCA = ggplot(df.PCA, aes(x = PC1, y = PC2, color = Species, shape = Technique, label = Biosample)) +
	geom_point(size = 2) +
	geom_text(size = 2) +
	theme_classic() +
	scale_color_manual(values = viridis(length(unique(df.PCA$Species))+1, option = "turbo")) +
	theme(axis.text.x = element_text(size =6),
    	  axis.text.y = element_text(size = 6),
    	  axis.title.x = element_text(size = 8),
    	  axis.title.y = element_text(size = 8),
    	  strip.text = element_text(size = 8),
          legend.position = "bottom",
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8))


P.AB = plot_grid(P.cor.bottom, P.PCA, nrow = 1, labels = c("A", "B"), rel_widths = c(1.3, 1), label_size = 12)


ggsave("new_plots/Figure_SX_Correllation_clustering_on_PDAL-Seq_reads_in_allignmen.svg", P.AB, width = 6.5, height = 4, units = "in")

ggsave("new_plots/Figure_SX_pairwise_correlation_by_sample.svg", P.cor.top, width = 6.5, height = 4)

