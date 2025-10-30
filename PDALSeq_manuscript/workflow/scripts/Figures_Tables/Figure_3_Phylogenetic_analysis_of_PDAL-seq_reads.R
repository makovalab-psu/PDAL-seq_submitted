
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
	      axis.text.x = element_text(angle = 90)) +
	scale_y_continuous(limits = c(-0.5,1), breaks = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), labels = c("-0.5", "", "0", "", "0.5", "", "1")) +
	scale_color_manual(values = c("red", "blue", "grey"))

ggsave("plots/Figure_SX_pairwise_correlation_by_sample/Figure_SX.svg", P.cor.top, width = 7, height = 6)

P.cor.bottom = ggplot(df.cor, aes(x = Sample_1, y = Sample_2, fill = spearman)) +
	geom_tile() +
	scale_fill_viridis(option = "turbo") +
	theme_classic() +
	theme(axis.title = element_blank(),
	      axis.text.x = element_blank()) 
	     

####B PCA plot####

df.PCA = read.csv("data/07_PhyloHGMP_model/alignment_clustering/hominoidae/PDALseq_PCA_in_alignment.csv")

df.PCA$Species = factor(df.PCA$Species, levels = unique(df.PCA$Species))

P.PCA = ggplot(df.PCA, aes(x = PC1, y = PC2, color = Species, shape = Technique, label = Biosample)) +
	geom_point(size = 4) +
	geom_text() +
	theme_classic() +
	scale_color_manual(values = viridis(length(unique(df.PCA$Species))+1, option = "turbo"))

####D-? Hidden Markov model#####

#Plot that shows signal in each state

df = readr::read_csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/Average_PDALseq_in_GHMM_states.csv") %>% arrange(MEAN)
State.levels = unique(df$state)

df$state = factor(df$state, levels = State.levels)
df$Experiment = factor(df$Experiment, levels = Exp.levels, labels = Exp.label)

P.tile = ggplot(df, aes(x = Experiment, y = state, fill = MEAN)) +
	geom_tile() +
	scale_fill_gradient(low = "white", high = "red") +
	theme_classic() +
	theme(axis.text.x = element_blank(),
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

df.count$Million_bp = round(df.count$count/10^6, digits = 0)

write.csv(df.count %>% arrange(state), "tables/Table_1_state_BP_counts.csv", row.names = F)

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

####make a plot that shows functional annotation enrichment in each state####

df.ann = read.csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_functional_element_enrichment/hsapiens.csv")

ANNOTATIONS = c("promoter_protCode_singleTrans",   "utr5_protCode_singleTrans", "gene_protCode_singleTrans", "gene_nonprotCode", "cds_protCode_singleTrans", "exon_protCode_singleTrans", "intron_protCode_singleTrans", "utr3_protCode_singleTrans", "enhancer_protCode_singleTrans", "rna_protCode_singleTrans", "chrG_cpgislands", "core_stochastic_origins_chm13v2.0_hglft.ucsc", "core_origins_hglft_genome_e0596_1cc9e0", "stochastic_origins_hglft_genome_2ff9ed_1ccfd0", "Recombination",   "repeats_protCode_singleTrans", "ngnr_protCode_singleTrans")

ANNOTATIONS.labels = c("promoter", "5'-UTR", "gene", "gene not protein coding", "cds", "exon", "intron", "3'-UTR",  "enhancer", "RNA", "CpG islands", "Replication origins", "Core origins", "Stochastic origins", "Recombination sites",   "Repeats" , "NGNR")

df.ann$State = factor(df.ann$State, levels = State.levels)
df.ann$ann = factor(df.ann$ann, levels = ANNOTATIONS, labels = ANNOTATIONS.labels)
df.ann$Sig = rep("ns", nrow(df.ann))
df.ann$Sig[which(df.ann$Wilcox.p <= 0.05)] = "*"
df.ann$Sig[which(df.ann$Wilcox.p <= 0.01)] = "**"

df.ann = df.ann %>% filter(ann %in% c("promoter", "enhancer",  "5'-UTR", "exon", "intron", "3'-UTR", "Repeats", "NGNR", "Recombination sites", "Core origins", "Stochastic origins"))
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

####Make a plot that show methylation in each state####

df.5mC = read.csv("results/10_Methylation_analysis/Methylation_by_state_and_species.csv")
df.5mC$State = factor(df.5mC$State, levels = rev(State.levels))

df.5mC$Species = factor(df.5mC$Species, levels = rev(c("Hsapien","Ptroglodytes", "Ppaniscus",    "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus")), labels = rev(c("H. sapiens","P. troglodytes", "P. paniscus",    "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus")))

P.5mC = ggplot(df.5mC, aes(x = Range, y = Species, fill = Percent)) +
	facet_wrap(~State, ncol = 1, strip.position = "left") +
	geom_tile() +
	scale_fill_gradient(low = "white", high = "black") +
	theme_classic() +
	scale_y_discrete(position = "right") +
	theme(axis.text.x = element_text(angle = 90),
	      legend.position = "right") +
	ggtitle("5mC rate at CpG sites") +
	xlab("5mC rate")
df.hypo = df.5mC %>% filter(Range == "0_0.2")	
df.hypo$State = factor(df.hypo$State, levels = State.levels)

library(ggbeeswarm)

P.hypo = ggplot(df.hypo , aes(x = Percent, y = State)) +
	geom_bar(stat="summary", fun = mean, color = "black", fill = "grey") +
    stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
    geom_beeswarm() +
    theme_classic() +
    ggtitle("Hypo-methylated\nCpG sites") +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank()) 
P.5mC = plot_grid(P.5mC, P.hypo, rel_widths = c(3,1), nrow = 1, align = "h")    
ggsave("plots/Figure_3_Phylogenetic_PDALseq/5mC_plot_5_windows.svg", P.5mC, width = 7, height = 5, scale = 1.2)
#ggsave("plots/Figure_3_Phylogenetic_PDALseq/5mC_plot_10_windows.svg", P.5mC, width = 7, height = 5, scale = 1.2)


####Look at signal around ORI####

df.OR = read.csv("results/11_Origin_analysis_analysis/PDALseq_signal_around_origins.csv")

df.OR$origin = factor(df.OR$origin, levels = c("core_origins_hglft_genome_e0596_1cc9e0", "stochastic_origins_hglft_genome_2ff9ed_1ccfd0"), labels = c("core origins\nn = 62892", "stochastic origins\nn = 254837"))
df.OR$state = factor(df.OR$state, levels = c("all", "4"), labels = c("all", "state 4"))

df.OR$sample = factor(df.OR$sample, levels = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq"), labels = c("GM24385", "GM24143", "HEK293T", "CCL86-Raji"))
df.OR$shuffle = factor(df.OR$shuffle, levels = c("unscrambled", "scrambled"))

df.OR = df.OR %>% filter(state != "state 4")

P.ORI = ggplot(df.OR, aes(x = kb_range, y = Average_PDALseq, linetype = shuffle, color = origin)) +
	geom_line() +
	theme_classic() +
	facet_wrap(~sample , scales = "free") +
	scale_color_manual(values = c("black", "red")) +
	xlab("Distance from ORI (kb)") +
	ylab("PDAL-seq reads")

####Look at TSS signal around ORI####

df.TSS = read.csv("results/12_TSS_analysis_analysis/PDALseq_signal_around_TSS.csv")

df.TSS$state = factor(df.TSS$state, levels = c("all", "4", "7", "3", "1", "6", "0", "5","2"))

df.TSS$sample = factor(df.TSS$sample, levels = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq"), labels = c("GM24385", "GM24143", "HEK293T", "CCL86-Raji"))
df.TSS$shuffle = factor(df.TSS$shuffle, levels = c("unscrambled", "scrambled"))


P.TSS = ggplot(df.TSS, aes(x = kb_range - 1, y = Average_PDALseq, linetype = shuffle, color = sample)) +
	geom_line() +
	theme_classic() +
	facet_wrap(~state, scales = "free") +
	scale_color_manual(values = viridis(5)) +
	xlab("Distance from TSS (kb)") +
	ylab("PDAL-seq reads/(Million reads total)") +
	xlim(-25, 25)

ggsave("plots/Figure_3_Phylogenetic_PDALseq/TSS.svg", P.TSS, height = 6, width = 14)

####Plot hmm-States Genome wide#### 

df.bed = readr::read_delim("resources/alignments/hominoidae/Hsapien.bed", delim = "\t", col_names = c("chr", "start", "end"))[,1:3]
df.state = readr::read_csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/GHMM_states.csv")
df.bed$hidden_state = df.state$hidden_state
df.bed$hidden_state = factor(df.bed$hidden_state, levels = c("2", "5", "0", "6", "1", "3", "7", "4"))

CHR = unique(df.bed$chr)
df.bed$chr = factor(df.bed$chr, levels = CHR)
df.chr.sizes = readr::read_delim("resources/genomes/Hsapien.chrom.sizes", col_names = c("chr", "end"))
df.chr.sizes$chr = factor(df.chr.sizes$chr, levels = CHR)
df.chr.sizes = df.chr.sizes %>% filter(!is.na(chr))

df.cen = readr::read_delim("resources/centromeres/Hsapien.bed", col_names = c("chr", "start", "end", "ann"))
df.cen$chr = factor(df.cen$chr, levels = CHR)
df.cen = df.cen %>% filter(!is.na(chr))

df.cband = readr::read_delim("resources/genomes/chm13v2.0_cytobands_allchrs.bed", col_names = c("chr", "start", "end", "P", "band"))
	
l.df.bed = df.bed %>% group_by(chr) %>% group_split()

l.df = {}

for (i in 1:length(l.df.bed)){
	df.i = data.frame(l.df.bed[[i]] %>% group_by(hidden_state) %>% summarize(BP = sum(end - start + 1)))
	Size = df.chr.sizes$end[which(df.chr.sizes$chr == l.df.bed[[i]]$chr[1])]
	df.i.nonsyntenic = data.frame("hidden_state" = c("nonsyntenic"), "BP" = Size - sum(df.i$BP))
	df.i = bind_rows(df.i, df.i.nonsyntenic)
	df.i$Percent = 100*df.i$BP/Size	
	df.i$chr = l.df.bed[[i]]$chr[1]
	l.df[[i]] = df.i
}

df.state.percent = bind_rows(l.df)
df.state.percent$chr = factor(df.state.percent$chr, levels = CHR)
df.state.percent$hidden_state = factor(df.state.percent$hidden_state, levels = c("nonsyntenic", "2", "5", "0", "6", "1", "3", "7", "4")) 

State_colors = c("grey", "grey", "grey", "grey", "grey", "#882D71", "#DB5829", "#F4A637")

Band_colors = c("Red", "#E2E2E2", "#383838", "#AAAAAA", "#8D8D8D", "#717171", "pink", "pink")
df.cband$chr = factor(df.cband$chr, levels = CHR)

P.ideogram = ggplot() +
	geom_rect(data = df.bed, mapping = aes(xmin = -0.8, xmax = 0.2, ymin = start, ymax = end, fill = hidden_state)) + #Color genome by HMM state
	scale_fill_manual(values = c(State_colors, Band_colors)) +
	geom_rect(data = df.cband, mapping = aes(xmin = -1.5, xmax = -1.0, ymin = start, ymax = end, fill = band)) +
	facet_wrap(~chr, nrow = 1) + 
	scale_y_reverse() + 
	theme_void() +
	theme(legend.position = "none",
	axis.title.x = element_blank(),
	axis.text.x = element_blank(),
	strip.text.y = element_blank() , 
  	strip.background = element_blank(),
	axis.ticks.x=element_blank(),
	axis.line.x=element_blank())
	
P.percent = ggplot(df.state.percent %>% filter(hidden_state %in% c("3", "7", "4")), aes(x = chr, y = Percent, fill = hidden_state)) +
	facet_wrap(~chr, nrow = 1, scales = "free_x") + 
	ggtitle("Ape state\n by hsa locus") +
	geom_bar(stat = "sum") +
	theme_classic() +
	scale_fill_manual(values = c("#882D71", "#DB5829", "#F4A637")) +
	guides(fill = guide_legend(nrow = 1)) +
	theme(legend.position = "top")
	
####Compile plots####

P.AB = plot_grid(P.cor.bottom, P.PCA, nrow = 1, labels = c("A", "B"), rel_widths = c(1.3, 1))
P.C = plot_grid(P.count, P.tile, P.con, P.nonB.E, P.ann.E, P.RNA, align = "h", nrow = 1, rel_widths = c(3,6,2,6,8,3), labels = c("C"))
P.D = plot_grid(P.percent, P.percent, ncol = 1,  labels = c("D",""), rel_heights = c(1.25,3), align = "h")
P.EF = plot_grid(P.5mC, P.ORI, nrow = 1, labels = c("E", "F"))
P = plot_grid(P.AB, P.C, P.D, P.EF, ncol = 1, rel_heights = c(1,1.2,1.6, 1)) 

ggsave("plots/Figure_3_Phylogenetic_PDALseq/Figure_3.svg", P, width = 8, height = 12, units = "in", scale = 2.5)
ggsave("plots/Figure_3_Phylogenetic_PDALseq/Figure_3_ideogram.png", P.ideogram, width = c(20/2.5), height = c(5.5/2.5), units = "in", scale = 2.5, dpi = 1200)
