
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
		 
Censat.levels = c("active_HOR", "inactive_HOR", "diverged_monomeric_HOR", "bSat", "HSat1", "HSat2", "HSat3", "Other_satellite", "rDNA", "Transition_region")
Censat.labels = c("active HOR", "inactive HOR", "div/mon HOR", "bSat", "HSat1", "HSat2", "HSat3", "Other satellite", "rDNA", "Transition region")
Censat_y =  c("active HOR", "inactive HOR", "div/mon HOR", "bSat", "HSat1", "HSat2", "HSat3", "Other satellite", "Transition region", "rDNA")
Censat.colors = c("#891641", "#CF6732", "#F8DB5B", "#DCAFD0", "#18988B", "#323467", "#77A1BA", "#768661", "#000000", "#959595")

nonB.levels = c("GQ", "SLS", "DR", "STR", "Z", "TRI", "MR", "CRU", "IR", "APR")
nonB.labels = c("G4", "SLS", "DR", "STR", "Z", "TRI", "MR", "CRU", "IR", "APR")

ANNOTATIONS = c("promoter_protCode_singleTrans",   "utr5_protCode_singleTrans", "gene_protCode_singleTrans", "gene_nonprotCode", "cds_protCode_singleTrans", "exon_protCode_singleTrans", "intron_protCode_singleTrans", "utr3_protCode_singleTrans", "enhancer_protCode_singleTrans", "rna_protCode_singleTrans", "chrG_cpgislands", "core_stochastic_origins_chm13v2.0_hglft.ucsc", "core_origins_hglft_genome_e0596_1cc9e0", "stochastic_origins_hglft_genome_2ff9ed_1ccfd0", "Recombination",   "repeats_protCode_singleTrans", "ngnr_protCode_singleTrans")

ANNOTATIONS.labels = c("promoter", "5'-UTR", "gene", "gene not protein coding", "cds", "exon", "intron", "3'-UTR",  "enhancer", "RNA", "CpG islands", "replication origins", "core origins", "stochastic origins", "recombination sites",   "repeats" , "NGNR")

####Redistribution of reads in functional elements####

df.read = readr::read_csv("results/21_add_K562_to_analysis/Read_enrichment.csv")

df.read$ann = factor(df.read$Ann, levels = ANNOTATIONS, labels = ANNOTATIONS.labels)
df.read$Sig = rep("ns", nrow(df.read))
df.read$Sig[which(df.read$Wilcox.p <= 0.05)] = "*"
df.read$Sig[which(df.read$Wilcox.p <= 0.01)] = "**"

df.read = df.read %>% filter(ann %in% c("promoter", "enhancer",  "5'-UTR", "exon", "intron", "3'-UTR", "repeats", "NGNR", "recombination sites", "core origins", "stochastic origins"))
df.read$LABEL = paste(round(df.read$Enrichment, digits = 1), df.read$Sig, sep = "\n") 
df.read$Exp = factor(df.read$Exp, levels = rev(c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq", "Hsapien_K562_PDALseq")), labels = rev(c("GM24385", "GM24143", "HEK293T", "CCL86-Raji", "K562")))

df.read$log2Enrichment[which(df.read$log2Enrichment > log2(3))] = log2(3)
df.read$log2Enrichment[which(df.read$log2Enrichment < -log2(3))] = -log2(3)

P.read = ggplot(df.read, aes(x = ann, y = Exp, fill = log2Enrichment, label = LABEL)) +
	geom_tile() +
	scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-log2(3), log2(3))) +
                         theme_classic() +
	geom_text(size = 2) + 
	ggtitle("PDAL-Seq reads in human functional elements") + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
	      axis.title.x = element_blank(),
	      legend.position = "top",
	      legend.title = element_text(size = 8, color = "black"),
	      legend.text = element_text(size = 6, color = "black"),
	      axis.title.y = element_blank(),
	      axis.text.y = element_text(size = 6, color = "black"),
	      legend.key.size = unit(0.1, "in"),
    	      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

ggsave("revised_plots/Figure_S16_PDAL-Seq_reads_in_human_functional_elements.svg", P.read , width = 4, height = 3.5, units = "in")
