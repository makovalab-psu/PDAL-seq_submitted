
library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)

####A-B Correllation plots####

Exp.levels   = c("Hsapien_GM24385_PDALseq",
		 "Hsapien_GM24143_PDALseq",
		 "Hsapien_HEK293T_PDALseq",
		 "Hsapien_CCL86-Raji_PDALseq")

Exp.label   = c("GM24385 PDAL-seq",
		 "GM24143 PDAL-seq",
		 "HEK293T PDAL-seq",
		 "CCL86-Raji")
State.levels = c("7", "0", "1", "5", "3", "2", "4", "6", "9", "8")

####D-? Hidden Markov model#####

#Plot that hows signal in each state

df = readr::read_csv("data/09_Human_HMM/04_Multivariate_guassian/10/Average_PDALseq_in_GHMM_states.csv") %>% arrange(MEAN)

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
	      #axis.text.y = element_blank(),
	      axis.title.y = element_blank()
	     ) +
	ggtitle("Mean PDAL-seq signal in state")

####Make a plot that shows the number of blocks in each state####

df.count = read.csv("data/09_Human_HMM/04_Multivariate_guassian/10/BP_in_state.csv", header = F)
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
	ggtitle("Percent BP") +
	xlab("percent (%)") 

####make a plot that shows non-B motif enrichment in each state####

df.nonB = read.csv("data/09_Human_HMM/04_Multivariate_guassian/10/State_nonB_enrichment/hsapiens.csv")

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

df.ann = read.csv("data/09_Human_HMM/04_Multivariate_guassian/10/State_functional_element_enrichment/hsapiens.csv")

ANNOTATIONS = c("promoter_protCode_singleTrans",   "utr5_protCode_singleTrans", "gene_protCode_singleTrans", "gene_nonprotCode", "cds_protCode_singleTrans", "exon_protCode_singleTrans", "intron_protCode_singleTrans", "utr3_protCode_singleTrans", "enhancer_protCode_singleTrans", "rna_protCode_singleTrans", "chrG_cpgislands", "core_stochastic_origins_chm13v2.0_hglft.ucsc", "core_origins_hglft_genome_e0596_1cc9e0", "stochastic_origins_hglft_genome_2ff9ed_1ccfd0", "Recombination",   "repeats_protCode_singleTrans", "ngnr_protCode_singleTrans")

ANNOTATIONS.labels = c("promoter", "5'-UTR", "gene", "gene not protein coding", "cds", "exon", "intron", "3'-UTR",  "enhancer", "RNA", "CpG islands", "Replication origins", "Core origins", "Stochastic origins", "Recombination sites",   "Repeats" , "NGNR")

df.ann$State = factor(df.ann$State, levels = State.levels)
df.ann$ann = factor(df.ann$ann, levels = ANNOTATIONS, labels = ANNOTATIONS.labels)
df.ann$Sig = rep("ns", nrow(df.ann))
df.ann$Sig[which(df.ann$Wilcox.p <= 0.05)] = "*"
df.ann$Sig[which(df.ann$Wilcox.p <= 0.01)] = "**"

df.ann = df.ann %>% filter(ann %in% c("promoter", "enhancer",  "5'-UTR", "exon", "intron", "3'-UTR", "Repeats", "NGNR", "Recombination sites",  "Core origins", "Stochastic origins"))
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

df.intersect = bind_rows(df.intersect %>% select(human_state, ape_state, Percent), data.frame("human_state" = factor(HUMAN_states, levels = State.levels), ape_state = NA, "Percent" = Percent))

P.intersect = ggplot(df.intersect, aes(human_state, x = Percent, fill = ape_state)) +
	geom_bar(stat = "identity") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1),
	      legend.position = "top",
	      legend.title = element_blank(),
	      axis.title.y = element_blank(),
	      axis.text.y = element_blank()) +
	scale_fill_manual(values = c("#386350", "#00C992", "#1964B0", "#7BB0DF", "#B6DBFF", "#882D71", "#DB5829", "#F4A637")) +
	ggtitle("Ape state\nintersect") +
	xlab("percent (%)") +
	scale_x_continuous(limits = c(0, 100))

####Human states mapped to chromosomes####

df.bed = readr::read_delim("data/09_Human_HMM/01_Remove_sex_chromosomes/Hsapien_GM24385_PDALseq.bg", delim = "\t", col_names = c("chr", "start", "end", "reads"))[,1:3]
df.state = readr::read_csv("data/09_Human_HMM/04_Multivariate_guassian/10/GHMM_states.csv")
df.bed$hidden_state = df.state$hidden_state
df.bed$hidden_state = factor(df.bed$hidden_state, levels = State.levels)

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

P.ideogram = ggplot() +
	geom_rect(data = df.bed, mapping = aes(xmin = -0.8, xmax = 0.2, ymin = start, ymax = end, fill = hidden_state)) + #Color genome by HMM state
	scale_fill_manual(values = c("grey", "grey", "#000000", "grey", "#B6DBFF", "#7BB0DF","grey", "#F4A637",  "#894B45", "Red")) +
	geom_rect(data = df.chr.sizes, mapping = aes(xmin = -1.5, xmax = -1.0, ymin = 0, ymax = end), fill = "grey", color = "black") + #draw chromosome
	geom_rect(data = df.cen, mapping = aes(xmin = -1.6, xmax = -0.9, ymin = start, ymax = end), fill = "dimgrey", color = "black") + #draw chromosome
	facet_wrap(~chr, nrow = 1) + 
	theme_void() +
	theme(legend.position = "none",
	axis.title.x = element_blank(),
	axis.text.x = element_blank(),
	strip.text.y = element_blank() , 
  	strip.background = element_blank(),
	axis.ticks.x=element_blank(),
	axis.line.x=element_blank()) +
	scale_y_reverse()
	
P.percent = ggplot(df.state.percent %>% filter(hidden_state %in% c("1", "3", "2", "6", "9", "8")), aes(x = chr, y = Percent, fill = hidden_state)) +
	facet_wrap(~chr, nrow = 1, scales = "free_x") + 
	ggtitle("Ape state\n by hsa locus") +
	geom_bar(stat = "sum") +
	theme_classic() +
	scale_fill_manual(values = c("#000000", "#B6DBFF", "#7BB0DF", "#F4A637", "#894B45", "Red")) +
	guides(fill = guide_legend(nrow = 1)) +
	theme(legend.position = "top")

####Redistribution of reads in functional elements####

df.read = read.csv("data/12_PDALseq/Read_enrichment.csv")

df.read$ann = factor(df.read$Ann, levels = ANNOTATIONS, labels = ANNOTATIONS.labels)
df.read$Sig = rep("ns", nrow(df.read))
df.read$Sig[which(df.read$Wilcox.p <= 0.05)] = "*"
df.read$Sig[which(df.read$Wilcox.p <= 0.01)] = "**"

df.read = df.read %>% filter(ann %in% c("promoter", "enhancer",  "5'-UTR", "exon", "intron", "3'-UTR", "Repeats", "NGNR", "Recombination sites", "Core origins", "Stochastic origins"))
df.read$LABEL = paste(round(df.read$Enrichment, digits = 1), df.read$Sig, sep = "\n") 
df.read$Exp = factor(df.read$Exp, levels = rev(c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq")), labels = rev(c("GM24385", "GM24143", "HEK293T", "CCL86-Raji")))

P.read = ggplot(df.read, aes(x = ann, y = Exp, fill = Enrichment, label = LABEL)) +
	geom_tile() +
	geom_text(color = "white") + 
	scale_fill_viridis(option = "turbo") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90)) +
	ggtitle("PDAL-seq reads in human functional elements") + 
	theme(axis.title.x = element_blank(),
	      legend.position = "top",
	      legend.title = element_blank(),
	      axis.title.y = element_blank())

####RM family enrichment in states####

#Intersect

df.states = read.csv("data/12_PDALseq/Human_state_rm_intersect_bp.csv")
df.baseline = read.csv("data/12_PDALseq/baseline_human_state_rm_intersect_bp.csv")
df.intersect = bind_rows(df.states, df.baseline)
df.intersect$Percent = 100*df.intersect$Intersect_bp/df.intersect$Total_bp_state


l.Class = strsplit(as.character(df.intersect$Ann), split = "/") 

Parent_class = c()
for(i in 1:nrow(df.intersect)){
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

df.intersect$Parent_class = factor(Parent_class, levels = Parent_class_levels)

df.intersect$State = factor(as.character(df.intersect$State), levels = c("Human HMM Region", paste("state_", State.levels, sep = "")))

P.RM.intersect = ggplot(df.intersect, aes(y = State, x = Percent)) +
	geom_bar(stat = "identity") +
	theme_classic() +
	geom_vline(data = df.intersect %>% filter(State == "Human HMM Region")  %>% summarize(SUM = sum(Percent)), mapping = aes(xintercept = SUM)) +
	#coord_flip() +
	theme(legend.position = "none") +
	ggtitle("Repeat content of\nhuman states")
#Enrichment
df.rm_E = read.csv("results/14_Repeat_masker_class_enrichment/03_RM_damily_enrichment.csv")

df.rm_E$State = factor(df.rm_E$State, levels = paste("state_", State.levels, sep = ""))
df.rm_E$ann = factor(df.rm_E$Ann, levels = c("Satelitte", "Low_complexity", "Simple_repeat", "Psuedogene", "DNA_transposons", "LTR", "SINE", "LINE"))
df.rm_E$Sig = rep("ns", nrow(df.rm_E))
df.rm_E$Sig[which(df.rm_E$Wilcox.p <= 0.05)] = "*"
df.rm_E$Sig[which(df.rm_E$Wilcox.p <= 0.01)] = "**"

df.rm_E$LABEL = paste(round(df.rm_E$Enrichment, digits = 1), df.rm_E$Sig, sep = "\n") 

P.rm_E = ggplot(df.rm_E, aes(x = ann, y = State, fill = Enrichment, label = LABEL)) +
	geom_tile() +
	geom_text(color = "white") + 
	scale_fill_viridis(option = "turbo") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90)) +
	ggtitle("Repeat masker annotation enrichment") + 
	theme(axis.title.x = element_blank(),
	      legend.position = "top",
	      legend.title = element_blank(),
	      axis.title.y = element_blank(),
	      axis.text.y = element_blank())

P.RM = plot_grid(P.RM.intersect, P.rm_E, nrow = 1, rel_widths = c(1,2), align = "h")
####Human state 8 GO analysis####

#Query
df.Hstate8.GO = read.delim("data/15_Gene_pathway_enrichments_Saswat/Human_state_8/goEnrichR_gene_GO.tsv", sep = "\t")
df.Hstate8.KEGG = read.delim("data/15_Gene_pathway_enrichments_Saswat/Human_state_8/goEnrichR_gene_KEGG.tsv", sep = "\t")
df.Hstate8.GO$Data_base = "GO"
df.Hstate8.KEGG$Data_base = "KEGG"
df.Hstate8 = bind_rows(df.Hstate8.GO, df.Hstate8.KEGG)

#Filter hits by significance and effect size####

df.Hstate8 = df.Hstate8 %>% filter(num_genes >= 25) %>% filter(adj_p_val <= 0.05)

#Null
df.Hstate9.GO = read.delim("data/15_Gene_pathway_enrichments_Saswat/Human_state_9/state9_randomised_GO_500times_with_genes.dat", sep = "\t")
df.Hstate9.KEGG = read.delim("data/15_Gene_pathway_enrichments_Saswat/Human_state_9/state9_randomised_KEGG_500times_with_genes.dat", sep = "\t")
df.Hstate9.GO$Data_base = "GO"
df.Hstate9.KEGG$Data_base = "KEGG"
df.Hstate9 = bind_rows(df.Hstate9.GO, df.Hstate9.KEGG)

#Test null that state 8 observations would be observed by randomly sampling state 9#

Obs.state9 = rep(NA, nrow(df.Hstate8))

for (i in 1:nrow(df.Hstate8)){
	Query = df.Hstate8$path_name[i]
	HITS = df.Hstate9$significant_occurences[which(df.Hstate9$X == Query)]
	if (length(HITS) == 0){
		Obs.state9[i] = 0
	}else{
		Obs.state9[i] = HITS
	}
}

df.Hstate8$Obs.state9 = Obs.state9
df.Hstate8$p.state9 = Obs.state9/100

#Write data for a supplemental table
df.Hstate8.table = df.Hstate8 %>% select(path_name, adj_p_val, num_genes, Data_base, Obs.state9, p.state9) 
write.csv(df.Hstate8.table, "tables/Table_SX_Human_state_8_Enrichr_analysis.csv", row.names = F)

#plot
df.Hstate8.plot = df.Hstate8.table %>% filter(p.state9 <= 0.05)
df.Hstate8.plot = df.Hstate8.plot[c(1:15, 40:50),]

####Short version for the human HMM figure####
df.Hstate8.plot$P.stat9.labels = paste("p =", as.character(df.Hstate8.plot$p.state9))
df.Hstate8.plot$P.stat9.labels[which(df.Hstate8.plot$P.stat9.labels == "p = 0")] = "p < 0.02"
df.Hstate8.plot$P.Enrichr.labels = paste("p =", as.character(round(df.Hstate8.plot$adj_p_val, digits = 2)))
df.Hstate8.plot$P.Enrichr.labels[which(df.Hstate8.plot$P.Enrichr.labels == "p = 0")] = "p < 0.02"


P.GO_out = ggplot(df.Hstate8.plot) +
	geom_rect(mapping = aes(xmax = num_genes, xmin = 0, ymin = (nrow(df.Hstate8.plot)-1):0, ymax = nrow(df.Hstate8.plot):1)) +
	geom_text(mapping = aes(x = -20, y = ((nrow(df.Hstate8.plot)-1):0 + nrow(df.Hstate8.plot):1)/2, label = P.stat9.labels)) +
	geom_text(mapping = aes(x = -60, y = ((nrow(df.Hstate8.plot)-1):0 + nrow(df.Hstate8.plot):1)/2, label = P.Enrichr.labels)) +
	theme_classic() +
	scale_x_continuous(position = "top", limits =c(-400,160), breaks = seq(-400, 150, by = 25)) +
	geom_text(mapping = aes(x = -250, y = ((nrow(df.Hstate8.plot)-1):0 + nrow(df.Hstate8.plot):1)/2, label = path_name)) +
	geom_text(mapping = aes(x = -400, y = ((nrow(df.Hstate8.plot)-1):0 + nrow(df.Hstate8.plot):1)/2, label = Data_base))
	

####Compile plots####

P.AB = plot_grid(P.count, P.tile, P.nonB.E, P.ann.E, P.intersect, align = "h", nrow = 1, rel_widths = c(2, 4, 6, 8, 8), labels = c("A", "", "", "", "B"))
P.ABC = plot_grid(P.AB, P.percent, ncol = 1, rel_heights = c(1,0.5), labels = c("", "C"))
P.DE = plot_grid(P.read, P.RM, nrow = 1, rel_widths = c(1,1.5), labels = c("D", "E")) 
P.ABCDE = plot_grid(P.ABC, P.DE, ncol = 1, rel_heights = c(1.75,1))
ggsave("plots/Figure_4_human_HMM_PDALseq/Figure_4ABCDE.svg", P.ABCDE, width = 8, height = 8, units = "in", scale = 2.5)
ggsave("plots/Figure_4_human_HMM_PDALseq/Figure_4C.png", P.ideogram, width = c(20/2.5), height = c(5.5/2.5), units = "in", scale = 2.5, dpi = 1200)
ggsave("plots/Figure_4_human_HMM_PDALseq/Figure_4D.svg", P.GO_out, width = 23/2.5, height = 10/2.5, units = "in", scale = 1.5)

system("touch plots/Figure_4_human_HMM_PDALseq/Figure_4.svg")

