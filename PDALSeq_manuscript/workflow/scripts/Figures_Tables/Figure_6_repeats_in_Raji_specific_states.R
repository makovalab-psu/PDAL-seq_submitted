
library(tidyverse)
library(viridis)
library(cowplot)

####Functional element enrichment####

ANNOTATIONS = c("promoter_protCode_singleTrans", "enhancer_protCode_singleTrans",  "utr5_protCode_singleTrans", "gene_protCode_singleTrans", "gene_nonprotCode", "cds_protCode_singleTrans", "exon_protCode_singleTrans", "intron_protCode_singleTrans", "utr3_protCode_singleTrans", "ngnr_protCode_singleTrans", "rna_protCode_singleTrans",   "repeats_protCode_singleTrans", "chrG_cpgislands", "Recombination", "core_stochastic_origins_chm13v2.0_hglft.ucsc", "core_origins_hglft_genome_e0596_1cc9e0", "stochastic_origins_hglft_genome_2ff9ed_1ccfd0")

ANNOTATIONS.labels = c("promoter", "enhancer",  "5'-UTR", "gene", "gene not protein coding", "cds", "exon", "intron", "3'-UTR", "NGNR", "RNA",   "Repeats", "CpG islands", "Recombination sites", "Replication origins", "Core origins", "Stochastic origins")

####Repeat masker enrichment in states####

State.levels = c("Human HMM Region", "state_7", "state_0", "state_1", "state_5", "state_3", "state_2", "state_4", "state_6", "state_9", "state_8")

df.RM = read.csv("data/12_PDALseq/Human_state_rm_enrichment.csv")

l.Class = strsplit(as.character(df.RM$Ann), split = "/") 

Parent_class = c()
for(i in 1:nrow(df.RM)){
	Parent_class[i] = l.Class[[i]][1]
}

#                      #Parent class                  #DFAM tree
Parent_class_levels = c("Unknown", 
			 "Low_complexity",
			 "Satellite", "Simple_repeat", #Tandem Repeat
			 "tRNA", "srpRNA", "snRNA", "scRNA", "rRNA", #Interspersed_Repeat ->  Pseudogene -> RNA 
			 "DNA?", "DNA", "RC",  #Interspersed_Repeat -> Transposable Element -> Class II DNA Transposition
			 "LTR", "Retroposon", "SINE", "LINE") #Interspersed_Repeat -> Transposable Element -> Class I Retrotransposition

Parent_class_labels = c("Unknown", 
			 "Low complexity",
			 "Tandem Repeat", "Tandem Repeat", #Tandem Repeat
			 "Interspersed Repeat\nPseudogene", "Interspersed Repeat\nPseudogene", "Interspersed Repeat\nPseudogene", "Interspersed Repeat\nPseudogene", "Interspersed Repeat\nPseudogene", #Interspersed_Repeat ->  Pseudogene -> RNA 
			 "Interspersed_Repeat\nTransposable Element\nClass II DNA Transposition", "Interspersed Repeat\nTransposable Element\nClass II DNA Transposition", "Interspersed Repeat\nTransposable Element\nClass II DNA Transposition",  #Interspersed_Repeat -> Transposable Element -> Class II DNA Transposition
			 "Interspersed_Repeat\nTransposable Element\nClass I Retrotransposition\nLTR", "Interspersed_Repeat\nTransposable Element\nClass I Retrotransposition\nRetroposon", "Interspersed_Repeat\nTransposable Element\nClass I Retrotransposition\nSINE", "Interspersed_Repeat\nTransposable Element\nClass I Retrotransposition\nLINE")	#Interspersed_Repeat -> Transposable Element -> Class I Retrotransposition		
df.RM$Parent_class = factor(Parent_class, levels = Parent_class_levels)

df.arrange = df.RM %>% filter(State == "state_2") %>% arrange(Enrichment) %>% arrange(Parent_class)

df.RM$Repeat_class = factor(df.RM$Ann, levels = df.arrange$Ann)
df.RM$Sig = "ns"
df.RM$Sig[which(df.RM$Wilcox.p <= 0.05)] = "*"
df.RM$Sig[which(df.RM$Wilcox.p <= 0.01)] = "**"

df.RM$State = factor(as.character(df.RM$State), levels = State.levels)
P.RM.E = ggplot(df.RM, aes(x = Repeat_class, y = Enrichment, fill = Parent_class)) +
	geom_bar(stat = "identity") +
	geom_text(mapping = aes(x = Repeat_class, y = Enrichment + 0.1, label = Sig)) + 
	theme_classic() +
	facet_wrap(~State, nrow = 1) +
	geom_hline(yintercept = 2) +
	coord_flip() 

ggsave("plots/Figure_SX_Repeat_analysis_in_Raji_nonB_state/Figure_SX.svg", P.RM.E, width = 7, height = 9, units = "in", scale = 2)

####Intersect####

df.states = read.csv("data/12_PDALseq/Human_state_rm_intersect_bp.csv")
df.baseline = read.csv("data/12_PDALseq/baseline_human_state_rm_intersect_bp.csv")
df.intersect = bind_rows(df.states, df.baseline)
df.intersect$Percent = 100*df.intersect$Intersect_bp/df.intersect$Total_bp_state

df.intersect$Repeat_class = factor(df.intersect$Ann, levels = df.arrange$Ann)

l.Class = strsplit(as.character(df.intersect$Ann), split = "/") 

Parent_class = c()
for(i in 1:nrow(df.intersect)){
	Parent_class[i] = l.Class[[i]][1]
}

df.intersect$Parent_class = factor(Parent_class, levels = Parent_class_levels)

df.intersect = df.intersect %>% arrange(Repeat_class) %>% arrange(State)


df.intersect$State = factor(as.character(df.intersect$State), levels = State.levels)

A = ggplot(df.intersect, aes(y = State, x = Percent)) +
	geom_bar(stat = "identity") +
	theme_classic() +
	geom_vline(data = df.intersect %>% filter(State == "Human HMM Region")  %>% summarize(SUM = sum(Percent)), mapping = aes(xintercept = SUM)) +
	#coord_flip() +
	theme(legend.position = "none") +
	ggtitle("Repeat content of human states")
	
B = ggplot(df.intersect %>% filter(Parent_class %in% c("Low_complexity", "Satellite", "Simple_repeat", "DNA", "LTR", "Retroposon", "SINE", "LINE")), aes(y = State, x = Percent, fill = Parent_class)) +
	geom_bar(stat = "identity") +
	theme_classic() +
	geom_vline(data = df.intersect %>% filter(State == "Human HMM Region") %>% group_by(Parent_class) %>% summarize(SUM = sum(Percent)) %>% filter(Parent_class %in% c("Low_complexity", "Satellite", "Simple_repeat", "DNA", "LTR", "Retroposon", "SINE", "LINE")), mapping = aes(xintercept = SUM)) +
	facet_wrap(~Parent_class, nrow = 1) +
	theme(legend.position = "none",
	      scale.y.title = element_blank())
	
P.BC = plot_grid(A, B, nrow = 1, labels = c("B", "C"), rel_widths = c(1,3), align = "h")
P.ABC = plot_grid(P.ann.E, P.BC, labels = c("A", "C"), ncol = 1)	
####Figure S5#####

ggsave("plots/Figure_5_Repeat_analysis_in_Raji_nonB_state/Figure_5.svg", P.ABC, width = 7, height = 6, units = "in", scale = 2)

