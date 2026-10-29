
library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)

####A-B Correllation plots####
Key_size = 0.15
Exp.levels   = c("Hsapien_GM24385_PDALseq",
		 "Hsapien_GM24143_PDALseq",
		 "Hsapien_HEK293T_PDALseq",
		 "Hsapien_CCL86-Raji_PDALseq")

Exp.label   = c("GM24385 PDAL-seq",
		 "GM24143 PDAL-seq",
		 "HEK293T PDAL-seq",
		 "CCL86-Raji")
State.levels = c("7", "0", "1", "5", "3", "2", "4", "6", "9", "8")


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
NONB_labels = c("all", "APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS",  "G4")

####D-? Hidden Markov model#####

#Plot that hows signal in each state

df = readr::read_csv("data/09_Human_HMM/04_Multivariate_guassian/10/Average_PDALseq_in_GHMM_states.csv") %>% arrange(MEAN)

df$state = factor(df$state, levels = State.levels)
df$Experiment = factor(df$Experiment, levels = Exp.levels, labels = Exp.label)

P.tile = ggplot(df, aes(x = Experiment, y = state, fill = MEAN)) +
	geom_tile() +
	scale_fill_gradient(low = "white", high = "red", name = "Mean PDAL-seq signal") +
	theme_classic()+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
	      axis.title.x = element_blank(),
	      legend.position = "top",
	      legend.title = element_text(size = 8, color = "black"),
	      legend.text = element_text(size = 6, color = "black"),
	      axis.title.y = element_blank(),
	      axis.text.y = element_blank(),
	      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
	      legend.key.size = unit(0.1, "in")) 

####Make a plot that shows the number of blocks in each state####

Human_HMM_color = c("#BEBEBE", "#BEBEBE", "#000000", "#BEBEBE", "1964B0", "1964B0", "#BEBEBE", "#F4A637", "#894B45", "#D11F5E")

df.count = read.csv("data/09_Human_HMM/04_Multivariate_guassian/10/BP_in_state.csv", header = F)
colnames(df.count) = c("state", "count")
df.count$percent = 100*df.count$count/sum(df.count$count)
df.count$state = factor(df.count$state, levels = State.levels)

P.count = ggplot(df.count, aes(y = state, x = percent, fill = state)) +
	geom_bar(stat = "identity") +
	theme_classic() +
	scale_fill_manual(values = Human_HMM_color) +
	theme(axis.text.x = element_text(angle = 90, size = 6, color = "black"),
	      axis.title.x = element_text(size = 8, color = "black"),
	      legend.position = "none",
	      legend.title = element_text(size = 8, color = "black"),
	      axis.title.y =  element_text(size = 8, color = "black"),
	      axis.text.y =  element_text(size = 6, color = "black"),
	      plot.title =  element_text(size = 8, color = "black"),
	      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
	      legend.key.size = unit(0.1, "in")) +
	xlab("Percent HMM\nregion in\nstate(%)") +
	ylab("Hidden Markov State")

####make a plot that shows non-B motif enrichment in each state####

df.nonB = read.csv("data/09_Human_HMM/04_Multivariate_guassian/10/State_nonB_enrichment/hsapiens.csv")

NONB.select = c("APR", "CRU", "TRI", "Z", "DR", "G4")
df.nonB$log2Enrichment[which(df.nonB$log2Enrichment > log2(3))] = log2(3)
df.nonB$log2Enrichment[which(df.nonB$log2Enrichment < -log2(3))] = -log2(3)
df.nonB$nonB = factor(df.nonB$nonB, levels = rev(NONB_levels), labels = rev(NONB_labels))
df.nonB$State = factor(df.nonB$State, levels = State.levels)
df.nonB$Sig = rep("ns", nrow(df.nonB))
df.nonB$Sig[which(df.nonB$Wilcox.p <= 0.05)] = "*"
df.nonB$Sig[which(df.nonB$Wilcox.p <= 0.01)] = "**"
df.nonB = df.nonB %>% filter(nonB %in% NONB.select)
MINMAX = max(abs(df.nonB$Enrichment))
df.nonB$LABEL = paste(round(df.nonB$Enrichment, digits = 1), df.nonB$Sig, sep = "\n") 

P.nonB.E = ggplot(df.nonB, aes(x = nonB, y = State, fill = log2Enrichment, label = LABEL)) +
	geom_tile() +
	geom_text(color = "black", size = 2) + 
	scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-log2(3), log2(3)),
                         name = "non-B DNA motif\nlog2Enrichment") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
	      axis.title.x = element_blank(),
	      legend.position = "top",
	      legend.title = element_text(size = 8, color = "black"),
	      legend.text = element_text(size = 6, color = "black"),
	      axis.title.y = element_blank(),
	      axis.text.y = element_blank(),
	      legend.key.size = unit(0.1, "in"),
    	      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

####make a plot that shows functional motif enrichment in each state####

df.ann = read.csv("data/09_Human_HMM/04_Multivariate_guassian/10/State_functional_element_enrichment/hsapiens.csv")

ANNOTATIONS = c("promoter_protCode_singleTrans",   "utr5_protCode_singleTrans", "gene_protCode_singleTrans", "gene_nonprotCode", "cds_protCode_singleTrans", "exon_protCode_singleTrans", "intron_protCode_singleTrans", "utr3_protCode_singleTrans", "enhancer_protCode_singleTrans", "rna_protCode_singleTrans", "chrG_cpgislands", "core_stochastic_origins_chm13v2.0_hglft.ucsc", "core_origins_hglft_genome_e0596_1cc9e0", "stochastic_origins_hglft_genome_2ff9ed_1ccfd0", "Recombination",   "repeats_protCode_singleTrans", "ngnr_protCode_singleTrans")

ANNOTATIONS.labels = c("promoter", "5'-UTR", "gene", "gene not protein coding", "cds", "exon", "intron", "3'-UTR",  "enhancer", "RNA", "CpG islands", "replication origins", "core origins", "stochastic origins", "recombination sites",   "repeats" , "NGNR")

df.ann$State = factor(df.ann$State, levels = State.levels)
df.ann$ann = factor(df.ann$ann, levels = ANNOTATIONS, labels = ANNOTATIONS.labels)
df.ann$log2Enrichment[which(df.ann$log2Enrichment > log2(3))] = log2(3)
df.ann$log2Enrichment[which(df.ann$log2Enrichment < -log2(3))] = -log2(3)
df.ann$State = factor(df.ann$State, levels = State.levels)
df.ann$Sig = rep("ns", nrow(df.ann))
df.ann$Sig[which(df.ann$Wilcox.p <= 0.05)] = "*"
df.ann$Sig[which(df.ann$Wilcox.p <= 0.01)] = "**"

df.ann = df.ann %>% filter(ann %in% c("promoter", "enhancer",  "5'-UTR", "exon", "intron", "3'-UTR", "repeats", "NGNR", "recombination sites", "core origins", "stochastic origins"))
df.ann$LABEL = paste(round(df.ann$Enrichment, digits = 1), df.ann$Sig, sep = "\n") 

P.ann.E = ggplot(df.ann, aes(x = ann, y = State, fill = log2Enrichment, label = LABEL)) +
	geom_tile() +
	geom_text(color = "black", size = 2) + 
	scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-log2(3), log2(3)),
                         name = "Functional annotation\nlog2Enrichment") +
	theme_classic() + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
	      axis.title.x = element_blank(),
	      legend.position = "top",
	      legend.title = element_text(size = 8, color = "black"),
	      legend.text = element_text(size = 6, color = "black"),
	      axis.title.y = element_blank(),
	      axis.text.y = element_blank(),
	      legend.key.size = unit(0.1, "in"),
    	      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))


####Make a plot that shows human HMM states intersecting with ape HMM states####

df.intersect = read.csv("data/09_Human_HMM/04_Multivariate_guassian/10/intersect_results.csv")
l.human_state = strsplit(df.intersect$query, split = "/")
l.ape_state = strsplit(df.intersect$database, split = "/")
human_state = c()
ape_state = c()

for (i in 1:nrow(df.intersect)){
	human_state[i] = gsub(".bed", "", gsub("state_", "", l.human_state[[i]][5]))
	ape_state[i] = gsub(".bed", "", gsub("state_", "", l.ape_state[[i]][6]))
}

df.intersect$human_state = factor(human_state, levels = State.levels)
df.intersect$ape_state = factor(ape_state, levels = c("2", "5", "0", "6", "1", "3", "7", "4"))
df.intersect$Percent = 100*df.intersect$intersect/df.intersect$total

HUMAN_states = unique(df.intersect$human_state)
Percent = c() 

for (i in 1:length(HUMAN_states)){
	df.i = df.intersect %>% filter(human_state == HUMAN_states[i])
	Percent[i] = 100 - sum(df.i$Percent)
}

df.intersect1 = bind_rows(df.intersect %>% select(human_state, ape_state, Percent), data.frame("human_state" = factor(HUMAN_states, levels = State.levels), ape_state = NA, "Percent" = Percent))

P.intersect = ggplot(df.intersect1, aes(human_state, x = Percent, fill = ape_state)) +
	geom_bar(stat = "identity") +
	theme_classic() +
	scale_fill_manual(name = "Ape state\nintersect", values = c("#386350", "#386350", "#386350", "#B6DBFF", "#B6DBFF", "#882D71", "#DB5829", "#F4A637")) +
	xlab("Percent human state\nin ape states (%)") +
	scale_x_continuous(limits = c(0, 100)) + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
	      axis.title.x = element_text(size = 8, color = "black"),
	      legend.position = "top",
	      legend.title = element_text(size = 8, color = "black"),
	      legend.text = element_text(size = 6, color = "black"),
	      axis.title.y = element_blank(),
	      axis.text.y = element_blank(),
	      legend.key.size = unit(0.1, "in"),
    	      plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "in"))

####Human states mapped to chromosomes####

df.bed = readr::read_delim("data/09_Human_HMM/01_Remove_sex_chromosomes/Hsapien_GM24385_PDALseq.bg", delim = "\t", col_names = c("chr", "start", "end", "reads"))[,1:3]
df.state = readr::read_csv("data/09_Human_HMM/04_Multivariate_guassian/10/GHMM_states.csv")
df.bed$hidden_state = df.state$hidden_state
df.bed$hidden_state = factor(df.bed$hidden_state, levels = State.levels)

df.cband = readr::read_delim("resources/genomes/chm13v2.0_cytobands_allchrs.bed", col_names = c("chr", "start", "end", "P", "band"))
	
CHR = unique(df.bed$chr)
df.bed$chr = factor(df.bed$chr, levels = CHR)
df.chr.sizes = readr::read_delim("resources/genomes/Hsapien.chrom.sizes", col_names = c("chr", "end"))
df.chr.sizes$chr = factor(df.chr.sizes$chr, levels = CHR)
df.chr.sizes = df.chr.sizes %>% filter(!is.na(chr))

df.cen = readr::read_delim("resources/centromeres/Hsapien.bed", col_names = c("chr", "start", "end", "ann"))
df.cen$chr = factor(df.cen$chr, levels = CHR)
df.cen = df.cen %>% filter(!is.na(chr))
	
l.df.bed = df.bed %>% group_by(chr) %>% group_split()

l.df = {}

for (i in 1:length(l.df.bed)){
	df.i = data.frame(l.df.bed[[i]] %>% group_by(hidden_state) %>% summarize(BP = sum(end - start + 1)))
	Size = df.chr.sizes$end[which(df.chr.sizes$chr == l.df.bed[[i]]$chr[1])]
	df.i.nonsyntenic = data.frame("hidden_state" = c("blacklist"), "BP" = Size - sum(df.i$BP))
	df.i = bind_rows(df.i, df.i.nonsyntenic)
	df.i$Percent = 100*df.i$BP/Size	
	df.i$chr = l.df.bed[[i]]$chr[1]
	l.df[[i]] = df.i
}

df.state.percent = bind_rows(l.df)
df.state.percent$chr = factor(df.state.percent$chr, levels = CHR)
df.state.percent$hidden_state = factor(df.state.percent$hidden_state, levels = c("blacklist", State.levels)) 

CHR = "chr19"

Band_colors = c("white", "#E2E2E2", "#383838", "#AAAAAA", "#8D8D8D", "#717171", "white", "white")
df.cband$chr = factor(df.cband$chr, levels = CHR)
P.ideogram = ggplot() +
	geom_rect(data = df.bed %>% filter(chr == CHR)%>% filter(hidden_state == 9), mapping = aes(xmin = -0.8, xmax = 0.2, ymin = start, ymax = end, fill = hidden_state)) + #Color genome by HMM state
	scale_fill_manual(values = "#D11F5E") +
	geom_rect(data = df.chr.sizes %>% filter(chr == CHR), mapping = aes(xmin = -1.5, xmax = -1.0, ymin = 0, ymax = end), fill = "grey", color = "black") + #draw chromosome
	geom_rect(data = df.cen %>% filter(chr == CHR), mapping = aes(xmin = -1.6, xmax = -0.9, ymin = start, ymax = end), fill = "dimgrey", color = "black") + #draw chromosome
	facet_wrap(~chr, nrow = 1) + 
	theme_void()+
	theme(legend.position = "none",
	axis.title.x = element_blank(),
	axis.text.x = element_blank(),
	strip.text.y = element_blank() , 
  	strip.background = element_blank(),
	axis.ticks.x=element_blank(),
	axis.line.x=element_blank(),
	panel.spacing = unit(0,'lines'))

ggsave("plots/Ideogram_for_Sallys_grand.png", P.ideogram, width = 0.8, height = 6.2, units = "in", dpi = 1200)


