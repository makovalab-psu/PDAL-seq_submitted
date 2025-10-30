####Figure 1D Correlation####

library(tidyverse)
library(viridis)
library(cowplot)
library(ggpubr)
library(ggbeeswarm)

v.margins = c(0.1,0.1,0.1,0.1)

#data/03_Correlation_and_black_list_unmappable_regions/04_genome_window_coverage_for_blacklist/High_MapQ/Hsapien_GM24385_TruSeq.bg

####99% Winsorization####

df.windows = read.csv("results/04_Genome_wide_comparison_to_nonB_motifs/nonB_coverage_versus_read.csv") %>% filter(nonB == "all")
df.windows$Experiment = paste(df.windows$Species, df.windows$Bio_sample, df.windows$Technique, sep = "_")

####99% Winsorization####

l.df = df.windows %>% group_by(Experiment) %>% group_split()
l.df2 ={}
for (i in 1:length(l.df)){
    l.df2[[i]] = data.frame(l.df[[i]])
    Read.Q99 = quantile(l.df[[i]]$reads, 0.99)
    Motif.Q99 = quantile(l.df[[i]]$nonb.count, 0.99)
    l.df[[i]]$reads[which(l.df[[i]]$reads > Read.Q99)] = as.integer(Read.Q99)
    l.df[[i]]$nonb.count[which(l.df[[i]]$nonb.count > Motif.Q99)] = as.integer(Motif.Q99)
    l.df[[i]] = data.frame(l.df[[i]])
    l.df2[[i]]$reads = l.df[[i]]$reads/Read.Q99
    l.df2[[i]]$nonb.count = l.df[[i]]$nonb.count/Motif.Q99
    l.df2[[i]]$interval = 1:nrow(l.df2[[i]]) 
}


df.windows = bind_rows(l.df2) %>% filter(Technique != "TruSeq")

species = "Hsapien"
cen = "resources/centromeres/Hsapien.bed"
v.levels = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq")
v.labels = c("GM24385",  "GM24143", "HEK293T", "CCL86-Raji")
v.heights = c(0.25, 0.25, 1, 3)
sat.file = "results/04_Genome_wide_comparison_to_nonB_motifs/hsat_2_3_annotations/Hsapien.bed"
sat.title = "hsat2,3"

df.reads = df.windows %>% filter(Species == species)
df.reads$Experiment = factor(df.reads$Experiment, levels = v.levels, labels = v.labels) 

#Plot reads v coordinate
P.reads = ggplot(df.reads, aes(x = interval, y = reads, fill = chr)) +
       facet_wrap(~Experiment, ncol = 1) +
	geom_bar(stat = "identity") +
	theme_classic() +
	scale_fill_manual(values = viridis(length(unique(df.reads$chr)), option = "turbo")) +
	scale_x_continuous(limits = c(min(df.reads$interval), max(df.reads$interval)), expand = c(0, 0))  +
  	theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank()) +
          ylab("reads")
	
#Plot non-B v coordinate
	
df.nonb = read.csv("results/04_Genome_wide_comparison_to_nonB_motifs/nonB_coverage_versus_read.csv") 
df.nonb$Experiment = paste(df.nonb$Species, df.nonb$Bio_sample, df.nonb$Technique, sep = "_")

df.nonb = df.nonb %>% filter(Experiment == "Hsapien_GM24385_PDALseq")

####99% Winsorization####

l.df = df.nonb %>% group_by(nonB) %>% group_split()
l.df2 ={}
for (i in 1:length(l.df)){
    l.df2[[i]] = data.frame(l.df[[i]])
    Read.Q99 = quantile(l.df[[i]]$reads, 0.99)
    Motif.Q99 = quantile(l.df[[i]]$nonb.count, 0.99)
    l.df[[i]]$reads[which(l.df[[i]]$reads > Read.Q99)] = as.integer(Read.Q99)
    l.df[[i]]$nonb.count[which(l.df[[i]]$nonb.count > Motif.Q99)] = as.integer(Motif.Q99)
    l.df[[i]] = data.frame(l.df[[i]])
    l.df2[[i]]$reads = l.df[[i]]$reads/Read.Q99
    l.df2[[i]]$nonb.count = l.df[[i]]$nonb.count/Motif.Q99
    l.df2[[i]]$interval = 1:nrow(l.df2[[i]]) 
}

df.nonb = bind_rows(l.df2)

NONB = c("GQ", "SLS", "Z", "TRI", "CRU", "APR", "all")

df.nonb = df.nonb %>% filter(nonB %in% NONB)
df.nonb$nonB = factor(df.nonb$nonB, levels = NONB)

P.motif = ggplot(df.nonb , aes(x = interval, y = nonB, fill = nonb.count)) +
        geom_tile() +
	theme_classic() +
	scale_fill_gradient(low = "white", high = "red") +
	scale_x_continuous(limits = c(min(df.nonb$interval), max(df.nonb$interval)), expand = c(0, 0)) +
    	theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
          ylab("non-B motifs")
	
	#Label in centromere coordinates
	df.cen = read.delim(cen, header = F)
	colnames(df.cen) = c("chr", "start", "end", "annotation")
	
	l.df.nonb = df.nonb %>% filter(nonB == "all") %>% group_by(chr) %>% group_split()
          
	l.df.centromere = {}

	for (i in 1:length(l.df.nonb)){
		df.cen.i = df.cen %>% filter(chr == l.df.nonb[[i]]$chr[i])
		Centromere = rep(NA, nrow(l.df.nonb[[i]]))
		for (j in 1:nrow(df.cen.i)){
			Min.start = which.min(abs(l.df.nonb[[i]]$start - df.cen.i$start[j]))
			Min.end = which.min(abs(l.df.nonb[[i]]$end - df.cen.i$end[j]))
			if (length(Min.start) > 0){Centromere[Min.start:Min.end] = 1}
		}
		l.df.centromere[[i]] = data.frame(l.df.nonb[[i]])
		l.df.centromere[[i]]$Centromere = Centromere
		
	}

	df.centromere = bind_rows(l.df.centromere)
	
	#Label sat coordinates
	df.sat = read.delim(sat.file, header = F)
	colnames(df.sat) = c("chr", "start", "end")
	
	l.df.sat = {}
	
	for (i in 1:length(l.df.nonb)){
		df.sat.i = df.sat %>% filter(chr == l.df.nonb[[i]]$chr[i])
		Sat = rep(FALSE, nrow(l.df.nonb[[i]]))
		for (j in 1:nrow(df.sat.i)){
			Min.start = which.min(abs(l.df.nonb[[i]]$start - df.sat.i$start[j]))
			Min.end = which.min(abs(l.df.nonb[[i]]$end - df.sat.i$end[j]))
			if (length(Min.start) > 0){Sat[Min.start:Min.end] = TRUE}
		}
		l.df.sat[[i]] = data.frame(l.df.nonb[[i]])
		l.df.sat[[i]]$Sat = Sat
		
	}
	
	df.Sat = bind_rows(l.df.sat)
	
	#Make an ideogram
	P.ideogram = ggplot(df.centromere, aes(x = interval, y = 1, fill = chr)) +
        geom_tile() +
        geom_point(mapping = aes(x = interval, y = Centromere), size = 5, color = "grey") + 
	theme_classic() +
	scale_fill_manual(values = viridis(length(unique(df.centromere$chr)), option = "turbo")) +
	scale_x_continuous(limits = c(min(df.centromere$interval), max(df.centromere$interval)), expand = c(0, 0)) +
    	theme_void() +
    	theme(legend.position = "none") +
    	ggtitle("Ideogram")
    	
    	#Maka a satalite plot
    	P.sat = ggplot(df.Sat, aes(x = interval, y = 1, fill = Sat)) +
        geom_tile() +
        theme_classic() +
	scale_x_continuous(limits = c(min(df.Sat$interval), max(df.Sat$interval)), expand = c(0, 0)) +
    	theme_void() +
    	scale_fill_manual(values = c("white", "black")) +
    	theme(legend.position = "none") +
    	ggtitle(sat.title)
    	
    	#Compile final plot
    	P = plot_grid(P.ideogram, P.sat, P.motif, P.reads, align = "v", ncol = 1, rel_heights = v.heights)
    	


#H. sapiens
P.Hsapiens = P	


P.A = P

####Hsat arrays####

T2T = c("Hsapien_GM24385", "Ptroglodytes_AG18354", "Ppaniscus_KB8711", "Ggorilla_KB3781", "Ppygmaeus_AG05252", "Pabelii_AG06213", "Ssyndactylus_JAMBI")
T2T.labels = c("H. sapiens\nGM24385", "P. troglodytes\nAG18354", "P. paniscus\nKB8711", "G. gorilla\nKB3781", "P. pygmaeus\nAG05252", "P. abelii\nAG06213", "S. syndactylus\n JAMBI")

df.sat.count = read.csv("results/06_reference_free_repeat_analysis/Short_tandem_repeats_by_technical_replicate.csv") %>% filter(Repeat == "AATGG") %>% filter(Bio_sample %in% T2T)

df.sat.count$Bio_sample = factor(df.sat.count$Bio_sample, levels = T2T, labels = T2T.labels) 

unique(paste(df.sat.count$Bio_sample, df.sat.count$Technique, sep = " "))

df.sat.count$Experiment = factor(paste(df.sat.count$Bio_sample, df.sat.count$Technique, sep = " "), levels = c("H. sapiens\nGM24385 PDAL-seq", "H. sapiens\nGM24385 TruSeq", "P. troglodytes\nAG18354 PDAL-seq", "P. troglodytes\nAG18354 TruSeq", "P. paniscus\nKB8711 PDAL-seq", "P. paniscus\nKB8711 TruSeq", "G. gorilla\nKB3781 PDAL-seq", "G. gorilla\nKB3781 TruSeq", "P. pygmaeus\nAG05252 PDAL-seq", "P. pygmaeus\nAG05252 TruSeq", "P. abelii\nAG06213 PDAL-seq", "P. abelii\nAG06213 TruSeq", "S. syndactylus\n JAMBI PDAL-seq", "S. syndactylus\n JAMBI TruSeq"), labels = c("H. sapiens PDAL-seq", "H. sapiens TruSeq", "P. troglodytes PDAL-seq", "P. troglodytes TruSeq", "P. paniscus PDAL-seq", "P. paniscus TruSeq", "G. gorilla PDAL-seq", "G. gorilla TruSeq", "P. pygmaeus PDAL-seq", "P. pygmaeus TruSeq", "P. abelii PDAL-seq", "P. abelii TruSeq", "S. syndactylus PDAL-seq", "S. syndactylus TruSeq"))

P.sat.count = ggplot(df.sat.count , aes(y = density.kbp.per.Mbp, x = Experiment, fill = Technique, group = Technique)) +
    geom_bar(stat="summary", fun = mean, color = "black") +
    stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
    geom_jitter() +
    theme_classic() +
    theme() +
    ylab("Sequencing rate (kbp/Million reads)") + 
    ggtitle("Reference free\nHsat2,3 monomers") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
    	  axis.title.x = element_blank()) +
    scale_fill_manual(values = c("white", "grey"))
    

df.lm = read.csv("results/05_nonB_motif_cluster_analysis/linear_models.csv") %>% filter(Fit == "Deduplicated parameters")

df.lm$Species = factor(df.lm$Species, levels = c("Hsapien","Ptroglodytes", "Ppaniscus",    "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus"), labels = c("H. sapiens","P. troglodytes", "P. paniscus",    "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus")) 
df.lm$X = factor(df.lm$X, levels = c("(Intercept)", "length.N", "GC.content", "APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS",  "GQ"))

P.lm = ggplot(data = df.lm, mapping = aes(x = X, y = Estimate)) +
	facet_wrap(~Technique, ncol = 1) +
	geom_bar(stat = "summary", fun = mean, fill = "grey") +
	stat_summary(fun.data = mean_sd, geom = "errorbar", color = "black") +
	geom_beeswarm(mapping = aes(color = Species)) +
	theme_classic() +
	scale_color_manual(values = viridis(length(unique(df.lm$Species)), option = "turbo")) +
	geom_hline(yintercept = 0) +
	theme(axis.title.x = element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1)) 

P.BC = plot_grid(P.sat.count, P.lm, ncol = 1, rel_heights = c(1,1.5), labels = c("B", "C"))

P.out = plot_grid(P.A, P.BC, nrow = 1, rel_widths = c(0.7, 0.3), labels = c("A"))

ggsave("plots/Poster_Figure_1.svg", P.out, width = 19, height = 9.5, units = "in", scale = 1.0, limitsize = FALSE)
	


