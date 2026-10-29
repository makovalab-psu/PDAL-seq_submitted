
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

NONB_levels = c("APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS",  "GQ", "all")
NONB_labels = c("APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS",  "G4", "all")
####Hidden Markov model#####

#Plot that shows signal in each state

df = readr::read_csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/Average_PDALseq_in_GHMM_states.csv") %>% arrange(MEAN)
as.character(unique(df$state))
State.levels = c("2", "5", "0", "6", "1", "3", "7", "4") 
State.labels = rev(c("8-origin", "7-recomb", "6-gene", "5-variable", "4-variable", "3-low", "2-low", "1-low"))

df$state = factor(df$state, levels = State.levels, labels = State.labels)
df$Experiment = factor(df$Experiment, levels = Exp.levels, labels = Exp.label)

df_mean = df %>% group_by(state) %>% summarize(MEAN = mean(MEAN))
P.tile = ggplot(df_mean, aes(x = "", y = state, fill = MEAN)) +
	geom_tile() +
	scale_fill_gradient(low = "white", high = "red", name = "Mean PDAL-seq signal") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
	      axis.title.x = element_blank(),
	      legend.position = "top",
	      legend.title = element_text(size = 8, color = "black"),
	      legend.text = element_text(size = 6, angle = 45, color = "black"),
	      axis.title.y = element_blank(),
	      legend.key.size = unit(0.1, "in"),
    	      axis.text.y = element_blank(),
	      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) 

####Make a plot that shows the number of blocks in each state####

df.count = read.csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/BP_in_state.csv", header = F)
colnames(df.count) = c("state", "count")
df.count$percent = 100*df.count$count/sum(df.count$count)
df.count$state = factor(df.count$state, levels = State.levels, labels = State.labels)

P.count = ggplot(df.count, aes(y = state, x = percent, fill = state)) +
	geom_bar(stat = "identity") +
	theme_classic() +
	scale_fill_manual(values = c("#386350", "#386350", "#386350", "#B6DBFF", "#B6DBFF", "#882D71", "#DB5829", "#F4A637")) +
	theme(axis.text.x = element_text(angle = 90, size = 6, color = "black"),
	      axis.title.x = element_text(size = 8, color = "black"),
	      legend.position = "none",
	      legend.title = element_text(size = 8, color = "black"),
	      axis.title.y =  element_text(size = 8, color = "black"),
	      axis.text.y =  element_text(size = 6, color = "black"),
	      plot.title =  element_text(size = 8, color = "black"),
	      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
	ggtitle("State size") +
	xlab("Percent\nalignment\nin state (%)") +
	ylab("Hidden Markov State")

####Make a plot that shows conservation in each state####

df.con = read.csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/mean_state_conservation.csv")

df.con$States = factor(df.con$States, levels = State.levels, labels = State.labels)

P.con = ggplot(df.con, aes(y = States, x = "", fill = Mean_conservation)) +
	geom_tile() +
	theme_classic() +
	theme(axis.text.x = element_blank(),
	      axis.title.x = element_text(size = 8, color = "black"),
	      legend.position = "top",
	      legend.title = element_text(size = 8, color = "black"),
	      legend.text = element_text(size = 6, color = "black", angle = 45, hjust = 1),
	      axis.title.y = element_blank(),
	      axis.text.y = element_blank(),
	      legend.key.size = unit(0.1, "in"),
    	      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
	xlab("percent (%)") +
	scale_fill_viridis("Conservation\nscore", trans = "log10") 
	
####make a plot that shows non-B motif enrichment in each state####

df.nonB = read.csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_nonB_enrichment/hsapiens.csv")

NONB.select = c("APR", "CRU", "TRI", "Z", "DR", "G4")
df.nonB$nonB = factor(df.nonB$nonB, levels = rev(NONB_levels), labels = rev(NONB_labels))
df.nonB$log2Enrichment[which(df.nonB$log2Enrichment > log2(3))] = log2(3)
df.nonB$log2Enrichment[which(df.nonB$log2Enrichment < -log2(3))] = -log2(3)
df.nonB$State = factor(df.nonB$State, levels = State.levels, labels = State.labels)
df.nonB$Sig = rep("ns", nrow(df.nonB))
df.nonB$Sig[which(df.nonB$Wilcox.p <= 0.05)] = "*"
df.nonB$Sig[which(df.nonB$Wilcox.p <= 0.01)] = "**"
df.nonB = df.nonB %>% filter(nonB %in% NONB.select)

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

####make a plot that shows functional annotation enrichment in each state####

df.ann = read.csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_functional_element_enrichment/hsapiens.csv")

ANNOTATIONS = c("promoter_protCode_singleTrans",   "utr5_protCode_singleTrans", "gene_protCode_singleTrans", "gene_nonprotCode", "cds_protCode_singleTrans", "exon_protCode_singleTrans", "intron_protCode_singleTrans", "utr3_protCode_singleTrans", "enhancer_protCode_singleTrans", "rna_protCode_singleTrans", "chrG_cpgislands", "core_stochastic_origins_chm13v2.0_hglft.ucsc", "core_origins_hglft_genome_e0596_1cc9e0", "stochastic_origins_hglft_genome_2ff9ed_1ccfd0", "Recombination",   "repeats_protCode_singleTrans", "ngnr_protCode_singleTrans")

ANNOTATIONS.labels = c("promoter", "5'-UTR", "gene", "gene not protein coding", "cds", "exon", "intron", "3'-UTR",  "enhancer", "RNA", "CpG islands", "replication origins", "core origins", "stochastic origins", "recombination sites",   "repeats" , "NGNR")

df.ann$State = factor(df.ann$State, levels = State.levels)
df.ann$ann = factor(df.ann$ann, levels = ANNOTATIONS, labels = ANNOTATIONS.labels)
df.ann$log2Enrichment[which(df.ann$log2Enrichment > log2(3))] = log2(3)
df.ann$log2Enrichment[which(df.ann$log2Enrichment < -log2(3))] = -log2(3)
df.ann$State = factor(df.ann$State, levels = State.levels, labels = State.labels)
df.ann$Sig = rep("ns", nrow(df.ann))
df.ann$Sig[which(df.ann$Wilcox.p <= 0.05)] = "*"
df.ann$Sig[which(df.ann$Wilcox.p <= 0.01)] = "**"

df.ann = df.ann %>% filter(ann %in% c("promoter", "enhancer",  "5'-UTR", "exon", "intron", "3'-UTR", "repeats", "NGNR", "recombination sites", "core origins", "stochastic origins"))
df.ann$LABEL = paste(round(df.ann$Enrichment, digits = 1), df.ann$Sig, sep = "\n") 

MINMAX = max(abs(df.ann$Enrichment))
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

####Make a plot that show RNA expression in each state####

df.RNA = read.csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_RNA_expression/hsapiens.csv")
df.RNA$State = factor(df.RNA$State, levels = State.levels, labels = State.labels)

P.RNA = ggplot(df.RNA, aes(y = State, x = Enrichment)) +
	geom_bar(stat = "identity") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, size = 6, color = "black"),
	      axis.title.x = element_text(size = 8, color = "black"),
	      legend.position = "top",
	      legend.title = element_text(size = 8, color = "black"),
	      axis.title.y =  element_blank(),
	      axis.text.y =  element_blank(),
	      plot.title =  element_text(size = 8, color = "black"),
	      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
	ggtitle("Transcription") +
	xlab("RNA-Seq reads\nenrichment") +
	geom_vline(xintercept = 1) 

####Make a plot that shows ATAC-Seq enrichment in each state####

df.ATAC = read.csv("results/18_ATACseq_analysis/ATACseq_read_enrichment.csv")
df.ATAC$State = factor(df.ATAC$State, levels = State.levels, labels = State.labels)

P.ATAC = ggplot(df.ATAC, aes(y = State, x = Enrichment)) +
	geom_bar(stat = "summary", fun = mean, fill = "grey") +
	stat_summary(fun.data = mean_sd, geom = "errorbar", color = "black") +
	geom_vline(xintercept = 1) +
	geom_point(mapping = aes(color = Bio_sample, shape = Bio_sample)) + 
	theme_classic() +
	scale_color_manual(values = c("dimgrey", "black")) +
	theme(axis.text.x = element_text(angle = 90, size = 6, color = "black"),
	      axis.title.x = element_text(size = 8, color = "black"),
	      legend.position = "top",
	      legend.title = element_text(size = 8, color = "black"),
	      axis.title.y =  element_blank(),
	      axis.text.y =  element_blank(),
	      plot.title =  element_text(size = 8, color = "black"),
	      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
	ggtitle("ATAC-Seq") +
	xlab("ATAC-Seq reads\nenrichment")  

####Plot hmm-States Genome wide#### 

df.bed = readr::read_delim("resources/alignments/hominoidae/Hsapien.bed", delim = "\t", col_names = c("chr", "start", "end"))[,1:3]
df.state = readr::read_csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/GHMM_states.csv")
df.bed$hidden_state = df.state$hidden_state
df.bed$hidden_state = factor(df.bed$hidden_state, levels = State.levels, labels = State.labels)

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
df.state.percent$hidden_state = factor(df.state.percent$hidden_state, levels = c("nonsyntenic", State.labels))

State_colors = c("grey", "#386350", "#386350", "#386350", "#B6DBFF", "#B6DBFF", "#882D71", "#DB5829", "#F4A637")


Band_colors = c("white", "#E2E2E2", "#383838", "#AAAAAA", "#8D8D8D", "#717171", "white", "white")
df.cband$chr = factor(df.cband$chr, levels = CHR)

P.ideogram = ggplot() +
	geom_rect(data = df.bed, mapping = aes(xmin = -0.8, xmax = 0.2, ymin = start, ymax = end, fill = hidden_state)) + #Color genome by HMM state
	scale_fill_manual(values = c("#386350", "#386350", "#386350", "#B6DBFF", "#B6DBFF", "#882D71", "#DB5829", "#F4A637", Band_colors)) +
	geom_rect(data = df.cband, mapping = aes(xmin = -1.5, xmax = -1.0, ymin = start, ymax = end, fill = band)) +
	geom_rect(data = df.cband %>% filter(band %in% c("gvar", "stalk")), mapping = aes(xmin = -1.5, xmax = -1.0, ymin = start, ymax = end), fill = "white", color = "black", linewidth = 0.1) +
	geom_rect(data = df.cband %>% filter(band == "acen"), mapping = aes(xmin = -1.6, xmax = -0.9, ymin = start, ymax = end), fill = "dimgrey", color = "dimgrey") +
	facet_wrap(~chr, nrow = 1) + 
	scale_y_reverse() + 
	theme_void() +
	theme(legend.position = "none",
	axis.title.x = element_blank(),
	axis.text.x = element_blank(),
	strip.text.y = element_blank() , 
  	strip.background = element_blank(),
	axis.ticks.x=element_blank(),
	axis.line.x=element_blank(),
	panel.spacing = unit(0,'lines'))

P.percent = ggplot(df.state.percent , aes(x = chr, y = Percent, fill = hidden_state)) +
	facet_wrap(~chr, nrow = 1, scales = "free_x") + 
	ggtitle("Distribution of ape PDAL-Seq HMM states in human chromosomes") +
	ylab("Percent (%)") +
	geom_bar(stat = "sum") +
	theme_classic() +
	scale_y_continuous(limits = c(0,100)) +
	scale_fill_manual(values = State_colors) +
	theme(axis.text.x = element_blank(),
	      axis.title.x = element_blank(),
	      legend.position = "top",
	      legend.title = element_text(size = 8, color = "black"),
	      axis.title.y =  element_text(size = 8, color = "black"),
	      axis.text.y =  element_text(size = 6, color = "black"),
	      plot.title =  element_text(size = 8, color = "black"),
	      strip.text = element_text(size = 6, color = "black"),
	      legend.key.size = unit(0.1, "in"),
    	      panel.spacing = unit(0,'lines')) 

####Look at signal around ORI####

df.OR = read.csv("results/21_add_K562_to_analysis/PDALseq_signal_around_origins.csv")
unique(df.OR$sample)

df.OR$origin = factor(df.OR$origin, levels = c("core_origins_hglft_genome_e0596_1cc9e0", "stochastic_origins_hglft_genome_2ff9ed_1ccfd0"), labels = c("core origins\nn = 62892", "stochastic origins\nn = 254837"))
df.OR$state = factor(df.OR$state, levels = c("all", "4"), labels = c("all", "state 4"))

df.OR$sample = factor(df.OR$sample, levels = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq", "Hsapien_K562_PDALseq"), labels = c("GM24385-HG002", "GM24143-HG004", "HEK293T", "CCL86-Raji", "K562"))
df.OR$shuffle = factor(df.OR$shuffle, levels = c("unscrambled", "scrambled"))

df.OR = df.OR %>% filter(state != "state 4")

P.ORI = ggplot(df.OR, aes(x = kb_range, y = Average_PDALseq/1000, linetype = shuffle, color = origin)) +
	geom_line() +
	theme_classic() +
	facet_wrap(~sample , scales = "free", nrow = 1) +
	scale_color_manual(values = viridis(3, option = "inferno")) +
	xlab("Distance from ORI (kb)") +
	ylab("PDAL-Seq reads (103)") +
	guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) + 
    	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
	      axis.title.x = element_text(size = 8),
	      legend.position = "right",
	      legend.title = element_text(size = 8),
	      legend.text =  element_text(size = 6),
	      axis.title.y =  element_text(size = 8),
	      axis.text.y =  element_text(size = 6),
	      plot.title =  element_text(size = 8)) 


####Compile plots####

P.A = plot_grid(P.count, P.tile, P.con, P.nonB.E, P.ann.E, P.RNA, P.ATAC, align = "h", nrow = 1, rel_widths = c(3,1.6,1.6,5.2,9.2,3,3), labels = c("A"), label_size = 12) 
P.A
P.ABC = plot_grid(P.A, P.percent, P.ORI, ncol = 1, rel_heights = c(2, 1, 1),  labels = c("","B","C"), label_size = 12) 

ggsave("revised_plots/Figure_3.svg", P.ABC, width = 6.5, height = 5, units = "in", scale = 1.5)
ggsave("revised_plots/Figure_3_ideogram.png", P.ideogram, width = 6.5, height = 2, units = "in", dpi = 1200)

