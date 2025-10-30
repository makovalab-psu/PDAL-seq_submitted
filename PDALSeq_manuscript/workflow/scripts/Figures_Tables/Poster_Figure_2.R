library(tidyverse)
library(cowplot)
library(viridis)

####C-G Hidden Markov model#####

df = readr::read_csv("results/07_PhyloHGMP_model/alignment_clustering/hominoidea/10/GHMM_states.csv")

####Make a tidy formatted data frame####

l.df = {}

for (i in 1:(ncol(df)-1)){
	df.i = data.frame("state" = df[[15]], "reads" = df[[i]], "Experiment" = colnames(df)[i])
	l.df[[i]] = df.i
}

df.tidy = bind_rows(l.df)

####Make a dataframe that summarizes the mean read count in each state####

l.df = {}

for (i in 1:(ncol(df)-1)){
	df.i = data.frame("state" = df[[15]], "reads" = df[[i]])
	l.df[[i]] = data.frame(df.i %>% group_by(state) %>% summarize(MEAN = mean(reads)))
	l.df[[i]]$Experiment = colnames(df)[i] 
}

df.out = bind_rows(l.df)

####Make a data frame to count and arrange the states####

df.arrange = df.tidy %>% group_by(state) %>% summarize(MEAN = mean(reads), COUNT = n()) %>% arrange(MEAN)

####Order the states by mean PDALseq reads####

df.arrange$state = factor(df.arrange$state, levels = df.arrange$state)
df.out$state = factor(df.out$state, levels = df.arrange$state)

####Make a plot that shows the number of blocks in each state####

P.count = ggplot(df.arrange, aes(y = state, x = COUNT)) +
	geom_bar(stat = "identity") +
	theme_classic() +
	theme(axis.title.x = element_blank(),
	      legend.position = "top") +
	ggtitle("Alignment\nblocks")

df.out$Experiment = factor(df.out$Experiment, levels = unique(df.out$Experiment))

P.tile = ggplot(df.out, aes(x = Experiment, y = state, fill = MEAN)) +
	geom_tile() +
	scale_fill_gradient(low = "white", high = "red") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90),
	      axis.title.x = element_blank(),
	      legend.position = "top") +
	ggtitle("Mean PDAL-seq signal in state")
	
####make a plot that shows non-B motif enrichment in each state####

df.nonB = read.csv("results/07_PhyloHGMP_model/alignment_clustering/hominoidea/10/State_nonB_enrichment/hsapiens.csv")

NONB = c("all", "APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS",  "GQ")

df.nonB$nonB = factor(df.nonB$nonB, levels = NONB)
df.nonB$State = factor(df.nonB$State, levels = df.arrange$state)
df.nonB$Sig = rep("ns", nrow(df.nonB))
df.nonB$Sig[which(df.nonB$Wilcox.p <= 0.05)] = "*"
df.nonB$Sig[which(df.nonB$Wilcox.p <= 0.01)] = "**"

P.nonB.E = ggplot(df.nonB, aes(x = nonB, y = State, fill = log2Enrichment, label = Sig)) +
	geom_tile() +
	geom_text() + 
	scale_fill_viridis(option = "turbo", limits = c(-2.2, 2.2)) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90))+
	ggtitle("non-B DNA motif enrichment") + 
	theme(axis.title.x = element_blank(),
	      legend.position = "top")

####make a plot that shows functional motif enrichment in each state####

df.ann = read.csv("results/07_PhyloHGMP_model/alignment_clustering/hominoidea/10/State_functional_element_enrichment/hsapiens.csv")

ANNOTATIONS = c("promoter_protCode_singleTrans", "enhancer_protCode_singleTrans",  "utr5_protCode_singleTrans", "gene_protCode_singleTrans", "gene_nonprotCode", "cds_protCode_singleTrans", "exon_protCode_singleTrans", "intron_protCode_singleTrans", "utr3_protCode_singleTrans", "rna_protCode_singleTrans",   "repeats_protCode_singleTrans", "chrG_cpgislands", "ngnr_protCode_singleTrans")

ANNOTATIONS.labels = c("promoter", "enhancer",  "5'-UTR", "gene", "gene not protein coding", "cds", "exon", "intron", "3'-UTR", "RNA",   "Repeats", "CpG islands", "NFNR")

df.ann$State = factor(df.ann$State, levels = df.arrange$state)
df.ann$ann = factor(df.ann$ann, levels = ANNOTATIONS, labels = ANNOTATIONS.labels)
df.ann$Sig = rep("ns", nrow(df.ann))
df.ann$Sig[which(df.ann$Wilcox.p <= 0.05)] = "*"
df.ann$Sig[which(df.ann$Wilcox.p <= 0.01)] = "**"

df.ann = df.ann %>% filter(ann %in% c("promoter", "enhancer",  "5'-UTR", "exon", "intron", "3'-UTR", "Repeats", "NFNR"))

P.ann.E = ggplot(df.ann, aes(x = ann, y = State, fill = log2Enrichment, label = Sig)) +
	geom_tile() +
	geom_text() + 
	scale_fill_viridis(option = "turbo", limits = c(-2, 2)) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90)) +
	ggtitle("Functional element enrichment") + 
	theme(axis.title.x = element_blank(),
	      legend.position = "top")


####Compile plots####

P = plot_grid(P.count, P.tile, P.nonB.E, P.ann.E, nrow = 1, align = "h", rel_widths = c(2,5,3,3))

ggsave("plots/Poster_Figure_2.svg", P, width = 20, height = 9.5, units = "in", scale = 1, limitsize = FALSE)
