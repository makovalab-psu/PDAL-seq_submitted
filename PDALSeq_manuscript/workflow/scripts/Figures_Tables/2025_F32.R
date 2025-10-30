
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
	theme(axis.title = element_blank(),
	      axis.text.x = element_blank()) +
	scale_y_continuous(limits = c(-0.5,1), breaks = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), labels = c("-0.5", "", "0", "", "0.5", "", "1")) +
	scale_color_manual(values = c("red", "blue", "grey"))

P.cor.bottom = ggplot(df.cor, aes(x = Sample_1, y = Sample_2, fill = spearman)) +
	geom_tile() +
	scale_fill_viridis(option = "turbo") +
	theme_classic() +
	theme(axis.title = element_blank(),
	      axis.text.x = element_blank()) 
	     
P.AB = plot_grid(P.cor.top, P.cor.bottom, ncol = 1, rel_heights = c(1,3), labels = c("A"), align = "v")

####B PCA plot####

df.PCA = read.csv("data/07_PhyloHGMP_model/alignment_clustering/hominoidae/PDALseq_PCA_in_alignment.csv")

df.PCA$Species = factor(df.PCA$Species, levels = unique(df.PCA$Species))

P.PCA = ggplot(df.PCA, aes(x = PC1, y = PC2, color = Species, shape = Technique, label = Biosample)) +
	geom_point(size = 4) +
	geom_text() +
	theme_classic() +
	scale_color_manual(values = viridis(length(unique(df.PCA$Species))+1, option = "turbo"))

####D-? Hidden Markov model#####

#Plot that hows signal in each state

df = readr::read_csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/Average_PDALseq_in_GHMM_states.csv") %>% arrange(MEAN)
State.levels = unique(df$state)

df$state = factor(df$state, levels = State.levels)
df$Experiment = factor(df$Experiment, levels = Exp.levels, labels = Exp.label)

P.tile = ggplot(df, aes(x = Experiment, y = state, fill = MEAN)) +
	geom_tile() +
	scale_fill_gradient(low = "white", high = "red") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90),
	      axis.title.x = element_blank(),
	      legend.position = "top",
	      legend.title = element_blank(),
	      axis.title.y = element_blank(),
	      axis.text.y = element_blank()) +
	ggtitle("Mean PDAL-seq signal in state")

####Make a plot that shows the number of blocks in each state####

df.count = read.csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/BP_in_state.csv", header = F)
colnames(df.count) = c("state", "count")
df.count$percent = df.count$count/sum(df.count$count)
df.count$state = factor(df.count$state, levels = State.levels)

P.count = ggplot(df.count, aes(y = state, x = percent)) +
	geom_bar(stat = "identity") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1),
	      axis.title.x = element_blank(),
	      legend.position = "top",
	      legend.title = element_blank()) +
	ggtitle("Percent BP in\nalignment") +
	xlab("percent (%)") 


####Make a plot that shows conservation in each state####

df.con = read.csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/mean_state_conservation.csv")

df.con$States = factor(df.con$States, levels = State.levels)

P.con = ggplot(df.con, aes(y = States, x = "", fill = Mean_conservation)) +
	geom_tile() +
	theme_classic() +
	theme(axis.text.x = element_blank(),
	      axis.title.x = element_blank(),
	      legend.position = "top",
	      legend.title = element_blank(),
	      axis.title.y = element_blank(),
	      axis.text.y = element_blank()) +
	ggtitle("Conservation") +
	xlab("percent (%)") +
	scale_fill_viridis(trans = "log10") 
	
####make a plot that shows non-B motif enrichment in each state####

df.nonB = read.csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_nonB_enrichment/hsapiens.csv")

NONB = c("all", "APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS",  "GQ")
NONB.select = c("APR", "CRU", "TRI", "Z", "DR", "GQ")
df.nonB$nonB = factor(df.nonB$nonB, levels = NONB)
df.nonB$State = factor(df.nonB$State, levels = State.levels)
df.nonB$Sig = rep("ns", nrow(df.nonB))
df.nonB$Sig[which(df.nonB$Wilcox.p <= 0.05)] = "*"
df.nonB$Sig[which(df.nonB$Wilcox.p <= 0.01)] = "**"
df.nonB = df.nonB %>% filter(nonB %in% NONB.select)
MINMAX = max(abs(df.nonB$Enrichment))
df.nonB$LABEL = paste(round(df.nonB$Enrichment, digits = 1), df.nonB$Sig, sep = "\n") 
P.nonB.E = ggplot(df.nonB, aes(x = nonB, y = State, fill = Enrichment, label = LABEL)) +
	geom_tile() +
	geom_text(color = "white") + 
	scale_fill_viridis(option = "turbo", limits = c(0, MINMAX)) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90))+
	ggtitle("non-B DNA motif enrichment") + 
	theme(axis.title.x = element_blank(),
	      legend.position = "top",
	      legend.title = element_blank(),
	      axis.title.y = element_blank(),
	      axis.text.y = element_blank())


####make a plot that shows functional motif enrichment in each state####

df.ann = read.csv("data/07_PhyloHGMP_model/alignment_clustering/hominoidae/8/State_functional_element_enrichmen/hsapiens.csv")

ANNOTATIONS = c("promoter_protCode_singleTrans", "enhancer_protCode_singleTrans",  "utr5_protCode_singleTrans", "gene_protCode_singleTrans", "gene_nonprotCode", "cds_protCode_singleTrans", "exon_protCode_singleTrans", "intron_protCode_singleTrans", "utr3_protCode_singleTrans", "ngnr_protCode_singleTrans", "rna_protCode_singleTrans",   "repeats_protCode_singleTrans", "chrG_cpgislands")

ANNOTATIONS.labels = c("promoter", "enhancer",  "5'-UTR", "gene", "gene not protein coding", "cds", "exon", "intron", "3'-UTR", "NGNR", "RNA",   "Repeats", "CpG islands")

df.ann$State = factor(df.ann$State, levels = State.levels)
df.ann$ann = factor(df.ann$ann, levels = ANNOTATIONS, labels = ANNOTATIONS.labels)
df.ann$Sig = rep("ns", nrow(df.ann))
df.ann$Sig[which(df.ann$Wilcox.p <= 0.05)] = "*"
df.ann$Sig[which(df.ann$Wilcox.p <= 0.01)] = "**"

df.ann = df.ann %>% filter(ann %in% c("promoter", "enhancer",  "5'-UTR", "exon", "intron", "3'-UTR", "Repeats", "NGNR"))
df.ann$LABEL = paste(round(df.ann$Enrichment, digits = 1), df.ann$Sig, sep = "\n") 

MINMAX = max(abs(df.ann$Enrichment))
P.ann.E = ggplot(df.ann, aes(x = ann, y = State, fill = Enrichment, label = LABEL)) +
	geom_tile() +
	geom_text(color = "white") + 
	scale_fill_viridis(option = "turbo", limits = c(0, MINMAX)) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90)) +
	ggtitle("Functional element enrichment") + 
	theme(axis.title.x = element_blank(),
	      legend.position = "top",
	      legend.title = element_blank(),
	      axis.title.y = element_blank(),
	      axis.text.y = element_blank())

####make a plot that shows recombination site enrichment in each state####

df.Recomb = read.csv("data/07_PhyloHGMP_model/alignment_clustering/hominoidae/8/State_functional_refseq_enrichment/hsapiens.csv")

Recomb = c("meiotic_recombination_region", "mitotic_recombination_region", "non_allelic_homologous_recombination_region", "recombination_feature")

Recomb.labels = c("meiotic", "mitotic", "non-allelic homologous", "recombination feature")

df.Recomb$State = factor(df.Recomb$State, levels = State.levels)
df.Recomb$ann = factor(df.Recomb$ann, levels = Recomb, labels = Recomb.labels)
df.Recomb$Sig = rep("ns", nrow(df.Recomb))
df.Recomb$Sig[which(df.Recomb$Wilcox.p <= 0.05)] = "*"
df.Recomb$Sig[which(df.Recomb$Wilcox.p <= 0.01)] = "**"
df.Recomb$LABEL = paste(round(df.Recomb$Enrichment, digits = 1), df.Recomb$Sig, sep = "\n") 

df.Recomb = df.Recomb %>% filter(ann %in% c("meiotic", "mitotic", "non-allelic homologous", "recombination feature"))

MINMAX = max(abs(df.Recomb$Enrichment), na.rm = T)
P.Recomb = ggplot(df.Recomb, aes(x = ann, y = State, fill = Enrichment, label = LABEL)) +
	geom_tile() +
	geom_text(color = "white") + 
	scale_fill_viridis(option = "turbo", limits = c(0, MINMAX)) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90)) +
	ggtitle("Recombination\nenrichment") + 
	theme(axis.title.x = element_blank(),
	      legend.position = "top",
	      legend.title = element_blank(),
	      axis.title.y = element_blank(),
	      axis.text.y = element_blank())

####Make a plot that show methylation in each state####

df.5mC = read.csv("data/07_PhyloHGMP_model/alignment_clustering/hominoidae/8/methylation/hsapiens.csv")
df.5mC$state = factor(df.5mC$state, levels = State.levels)

P.5mC = ggplot(df.5mC, aes(x = as.character(methylation_rate), y = state, fill = percent)) +
	geom_tile() +
	scale_fill_gradient(low = "white", high = "black") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90),
	      legend.position = "top",
	      legend.title = element_blank(),
	      axis.title.y = element_blank(),
	      axis.text.y = element_blank()) +
	ggtitle("5mC rate at CpG sites") +
	xlab("5mC")

####Make a plot that show RNA expression in each state####

df.RNA = read.csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_RNA_expression/hsapiens.csv")
df.RNA$State = factor(df.RNA$State, levels = State.levels)
df.RNA$State = factor(df.RNA$State, levels = State.levels)
P.RNA = ggplot(df.RNA, aes(y = State, x = Enrichment)) +
	geom_bar(stat = "identity") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1),
	      axis.title.x = element_blank(),
	      legend.position = "top",
	      legend.title = element_blank(),
	      axis.title.y = element_blank(),
	      axis.text.y = element_blank()) +
	ggtitle("RNA-seq\nEnrichment") +
	xlab("percent (%)") +
	geom_vline(xintercept = 1)




####Compile plots####

P.ABC = plot_grid(P.AB, P.PCA, nrow = 1, labels = c("", "C"), rel_widths = c(1.3, 1))
P.D = plot_grid(P.count, P.tile, P.con, P.nonB.E, P.ann.E, P.Recomb, P.5mC, P.RNA, align = "h", nrow = 1, rel_widths = c(3,12,3,6,8,4,5,3), labels = c("","D"))
P = plot_grid(P.ABC, P.D, ncol = 1, rel_heights = c(1,1.4)) 


ggsave("plots/Figure_3_Phylogenetic_PDALseq/Figure_3.svg", P, width = 8, height = 6, units = "in", scale = 2.5)

