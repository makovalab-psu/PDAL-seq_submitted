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

Male   = c("Hsapien_GM24385_PDALseq",
		 "Ptroglodytes_AG18354_PDALseq",
		 "Ppaniscus_KB8711_PDALseq",
		 "Ggorilla_KB3781_PDALseq",
		 "Ppygmaeus_AG05252_PDALseq",
		 "Pabelii_AG06213_PDALseq",
		 "Ssyndactylus_JAMBI_PDALseq")

TruSeq = c("Hsapien_GM24385_TruSeq",
		 "Ptroglodytes_AG18354_TruSeq",
		 "Ppaniscus_KB8711_TruSeq",
		 "Ggorilla_KB3781_TruSeq",
		 "Ppygmaeus_AG05252_TruSeq",
		 "Pabelii_AG06213_TruSeq",
		 "Ssyndactylus_JAMBI_TruSeq")
		 
Transformed = c("Hsapien_GM24385_PDALseq",
		 "Hsapien_GM24143_PDALseq",
		 "Ptroglodytes_AG18354_PDALseq",
		 "Ptroglodytes_AG18358_PDALseq",
		 "Ppaniscus_KB14571_PDALseq",
		 "Ggorilla_F7549_PDALseq",
		 "Pabelii_AG06213_PDALseq",
		 "Ssyndactylus_JAMBI_PDALseq",
		 "Ssyndactylus_Karenina_PDALseq")

Lymphoblastoid   = c("Hsapien_GM24385_PDALseq",
		 "Hsapien_GM24143_PDALseq",
		 "Hsapien_HEK293T_PDALseq",
		 "Hsapien_CCL86-Raji_PDALseq",
		 "Ptroglodytes_AG18354_PDALseq",
		 "Ptroglodytes_AG18358_PDALseq",
		 "Ssyndactylus_JAMBI_PDALseq",
		 "Ssyndactylus_Karenina_PDALseq")
		 
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

#Remove Raji/HEK293T, matching samples, and TruSeq

df = df.cor %>% filter(!Sample_1 %in% c("Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq")) %>% filter(!Sample_2 %in% c("Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq")) %>% filter(!Sample_1 == Sample_2)  %>% filter(!Sample_1 %in% TruSeq) %>% filter(!Sample_2 %in% TruSeq) 

####Deduplicate the distribution#####

PDALseq = c("Hsapien_GM24385_PDALseq",
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
		 "Ssyndactylus_Karenina_PDALseq")

Unique_comparisons = c()
for (i in 1:(length(PDALseq)-1)){
	for (j in (i+1):length(PDALseq)){
		Unique_comparisons = c(Unique_comparisons, paste(PDALseq[i], PDALseq[j]))
	}
}

nrow(df)
nrow(df[which(paste(df$Sample_1, df$Sample_2) %in% Unique_comparisons),])
df = df[which(paste(df$Sample_1, df$Sample_2) %in% Unique_comparisons),]

####Compare within species correlation to between species comparison####

l.Sample_1 = strsplit(df$Sample_1, split = "_")
l.Sample_2 = strsplit(df$Sample_2, split = "_")

Species_1 = c()
Species_2 = c()

for (i in 1:nrow(df)){
	Species_1[i] = l.Sample_1[[i]][1]
	Species_2[i] = l.Sample_2[[i]][1]
}

Species_levels = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus")
Species_labels = c("human", "chimp", "bonobo", "gorilla", "orangutan", "orangutan", "siamang")
Species_1 = factor(Species_1, levels = Species_levels, labels = Species_labels)
Species_2 = factor(Species_2, levels = Species_levels, labels = Species_labels)

Within_species = rep("Outside species", nrow(df))
Within_species[which(Species_1 == Species_2)] = "Inside species"


df$Within_species = Within_species



library(ggbeeswarm)

####Compare transformed to not transformed correllations####

Sample_1_transformed = rep("Primary", nrow(df))
Sample_1_transformed[which(df$Sample_1 %in% Transformed)] = "transformed"
Sample_2_transformed = rep("Primary", nrow(df))
Sample_2_transformed[which(df$Sample_2 %in% Transformed)] = "transformed"

Transformed_comparison = paste(Sample_1_transformed, Sample_2_transformed, sep = " to ")

Transformed_comparison = factor(Transformed_comparison, levels = c("Primary to Primary", "Primary to transformed", "transformed to Primary", "transformed to transformed"), labels = c("Primary to Primary", "Primary to transformed", "Primary to transformed", "Transformed to transformed"))

df$Transformed_comparison = Transformed_comparison

####Compare Lymphoblast to fibroblast correllations####

Sample_1_celltype = rep("Fibroblast", nrow(df))
Sample_2_celltype = rep("Fibroblast", nrow(df))
Sample_1_celltype[which(df$Sample_1 %in% Lymphoblastoid)] = "Lymphoblastoid"
Sample_2_celltype[which(df$Sample_2 %in% Lymphoblastoid)] = "Lymphoblastoid"

Cell_type_comparison = paste(Sample_1_celltype, Sample_2_celltype, sep = " to ")

Cell_type_comparison = factor(Cell_type_comparison, levels = c("Fibroblast to Fibroblast", "Fibroblast to Lymphoblastoid", "Lymphoblastoid to Fibroblast", "Lymphoblastoid to Lymphoblastoid"), labels = c("Fibroblast to Fibroblast", "Lymphoblastoid to Fibroblast", "Lymphoblastoid to Fibroblast", "Lymphoblastoid to Lymphoblastoid"))

df$Cell_type_comparison = Cell_type_comparison


####Compare sex correlations####

Sample_1_Sex = rep("Female", nrow(df))
Sample_2_Sex = rep("Female", nrow(df))
Sample_1_Sex[which(df$Sample_1 %in% Male)] = "Male"
Sample_2_Sex[which(df$Sample_2 %in% Male)] = "Male"

Sex_comparison = paste(Sample_1_Sex, Sample_2_Sex, sep = " to ")

Sex_comparison = factor(Sex_comparison, levels = c("Female to Female", "Female to Male", "Male to Female", "Male to Male"), labels = c("Female to Female", "Female to Male", "Female to Male", "Male to Male"))

df$Sex_comparison = Sex_comparison

####Make plots####

P.A = ggplot(df, aes(x = Within_species, y = spearman)) +
	stat_compare_means(method = "kruskal.test") +
	geom_bar(stat="summary", fun = mean, color = "black") +
    	stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
    	geom_beeswarm() +
	theme_classic() +
	scale_y_continuous(limits = c(0,1)) +
	theme(axis.title.x = element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1)) +
	ylab("Spearman's rank correlation") +
	ggtitle("Correlation by\nspecies")

P.B = ggplot(df, aes(x = Transformed_comparison, y = spearman)) +
	stat_compare_means(method = "kruskal.test") +
	geom_bar(stat="summary", fun = mean, color = "black") +
    	stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
    	geom_beeswarm() +
	theme_classic() +
	scale_y_continuous(limits = c(0,1)) +
	theme(axis.title.x = element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1)) +
	ylab("Spearman's rank correlation") +
	ggtitle("Correlation by\ncell type (Transformed/\nprimary)")

P.C = ggplot(df, aes(x = Cell_type_comparison, y = spearman)) +
	stat_compare_means(method = "kruskal.test") +
	geom_bar(stat="summary", fun = mean, color = "black") +
    	stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
    	geom_beeswarm() +
	theme_classic()+
	scale_y_continuous(limits = c(0,1)) +
	theme(axis.title.x = element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1)) +
	ylab("Spearman's rank correlation")+
	ggtitle("Correlation by\ncell type (lymphoblast/\nfibroblast)")

P.D = ggplot(df, aes(x = Sex_comparison, y = spearman)) +
	stat_compare_means(method = "kruskal.test") +
	geom_bar(stat="summary", fun = mean, color = "black") +
    	stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
    	geom_beeswarm() +
	theme_classic() +
	scale_y_continuous(limits = c(0,1)) +
	theme(axis.title.x = element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1)) +
	ylab("Spearman's rank correlation") +
	ggtitle("Correlation by\nsex")
	
P = plot_grid(P.A, P.B, P.C, P.D, nrow = 1, align = "h")
P
ggsave("plots/Figure_SX_correlation_by_species_cell_type_sex/Figure_SX.svg", P, width = 7, height = 4, units = "in", scale = 1.7)
