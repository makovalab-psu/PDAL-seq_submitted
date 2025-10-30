####Figure 1D Correlation####

library(tidyverse)
library(viridis)
library(cowplot)

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

df.windows = bind_rows(l.df2)

####Function to make plots by species####

Plot_reads_v_nonb = function(species = "Hsapien",
			      cen = "resources/centromeres/Hsapien.bed",
			      v.levels = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_CCL86-Raji_PDALseq", "Hsapien_HEK293T_PDALseq", "Hsapien_GM24385_TruSeq"),
			      v.labels = c("GM24385 PDAL-seq", "GM24143 PDAL-seq", "CCL86-Raji PDAL-seq", "HEK293T PDAL-seq", "GM24385 TruSeq"),
			      v.heights = c(0.25, 0.25, 1.2, 5),
			      sat.file = "results/04_Genome_wide_comparison_to_nonB_motifs/hsat_2_3_annotations/Hsapien.bed",
			      sat.title = "hsat2,3"){
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
          	ylab("PDAL-seq reads")
	
	#Plot non-B v coordinate
	
	df.nonb = read.csv("results/04_Genome_wide_comparison_to_nonB_motifs/nonB_coverage_versus_read.csv") %>% filter(Species == species)
	df.nonb$Experiment = paste(df.nonb$Species, df.nonb$Bio_sample, df.nonb$Technique, sep = "_")

	df.nonb = df.nonb %>% filter(Experiment == unique(df.nonb$Experiment)[1])
	
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
	
	NONB = c("GQ", "SLS", "DR", "STR", "Z", "TRI", "MR", "CRU", "IR", "APR", "all")
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
    	return(P)
} 


#H. sapiens
P.Hsapiens = Plot_reads_v_nonb()	

ggsave("plots/Figure_S2_Genome_wide_plots/Hsapien.svg", P.Hsapiens, width = 8, height = 5, units = "in", scale = 2.5)

#P. troglodytes

df = df.windows %>% filter(Species == "Ptroglodytes")

P.Ptroglodytes = Plot_reads_v_nonb(species = "Ptroglodytes",
			      cen = "resources/centromeres/Ptroglodytes.bed",
			      v.levels = c("Ptroglodytes_AG18354_PDALseq", "Ptroglodytes_AG18358_PDALseq", "Ptroglodytes_AG18354_TruSeq"),
			      v.labels = c("AG18354 PDAL-seq", "AG18358 PDAL-seq", "AG18354 TruSeq"),
			      v.heights = c(0.25, 0.25, 1.2, 3),
			      sat.file = "results/04_Genome_wide_comparison_to_nonB_motifs/hsat_2_3_annotations/Ptroglodytes.bed",
			      sat.title = "hsat2,3")	

ggsave("plots/Figure_S2_Genome_wide_plots/Ptroglodytes.svg", P.Ptroglodytes, width = 8, height = 3, units = "in", scale = 2.5)

#P. paniscus

df = df.windows %>% filter(Species == "Ppaniscus")

P.Ppaniscus = Plot_reads_v_nonb(species = "Ppaniscus",
			      cen = "resources/centromeres/Ppaniscus.bed",
			      v.levels = c("Ppaniscus_KB8711_PDALseq", "Ppaniscus_KB14571_PDALseq", "Ppaniscus_KB8711_TruSeq"),
			      v.labels = c("KB8711 PDAL-seq", "KB14571 PDAL-seq", "KB8711 TruSeq"),
			      v.heights = c(0.25, 0.25, 1.2, 3),
			      sat.file = "results/04_Genome_wide_comparison_to_nonB_motifs/hsat_2_3_annotations/Ppaniscus.bed",
			      sat.title = "hsat2,3")	

ggsave("plots/Figure_S2_Genome_wide_plots/Ppaniscus.svg", P.Ppaniscus, width = 8, height = 3, units = "in", scale = 2.5)
          

#G. gorilla

df = df.windows %>% filter(Species == "Ggorilla")

P.Ggorilla = Plot_reads_v_nonb(species = "Ggorilla",
			      cen = "resources/centromeres/Ggorilla.bed",
			      v.levels = c("Ggorilla_KB3781_PDALseq", "Ggorilla_F7549_PDALseq", "Ggorilla_KB3781_TruSeq"),
			      v.labels = c("KB3781 PDAL-seq", "F7549 PDAL-seq", "KB3781 TruSeq"),
			      v.heights = c(0.25, 0.25, 1.2, 3),
			      sat.file = "results/04_Genome_wide_comparison_to_nonB_motifs/hsat_2_3_annotations/Ggorilla.bed",
			      sat.title = "hsat2,3")	

ggsave("plots/Figure_S2_Genome_wide_plots/Ggorilla.svg", P.Ggorilla, width = 8, height = 3, units = "in", scale = 2.5)


#P. pygmaeus

df = df.windows %>% filter(Species == "Ppygmaeus")

P.Ppygmaeus = Plot_reads_v_nonb(species = "Ppygmaeus",
			      cen = "resources/centromeres/Ppygmaeus.bed",
			      v.levels = c("Ppygmaeus_AG05252_PDALseq", "Ppygmaeus_AG05252_TruSeq"),
			      v.labels = c("AG05252 PDAL-seq", "AG05252 TruSeq"),
			      v.heights = c(0.25, 0.25, 1.2, 2),
			      sat.file = "results/04_Genome_wide_comparison_to_nonB_motifs/hsat_2_3_annotations/Ppygmaeus.bed",
			      sat.title = "hsat2,3")	

ggsave("plots/Figure_S2_Genome_wide_plots/Ppygmaeus.svg", P.Ppygmaeus, width = 8, height = 2.5, units = "in", scale = 2.5)
	

#P. abelii

df = df.windows %>% filter(Species == "Pabelii")

P.Pabelii = Plot_reads_v_nonb(species = "Pabelii",
			      cen = "resources/centromeres/Pabelii.bed",
			      v.levels = c("Pabelii_AG06213_PDALseq", "Pabelii_AG06213_TruSeq"),
			      v.labels = c("AG06213 PDAL-seq", "AG06213 TruSeq"),
			      v.heights = c(0.25, 0.25, 1.2, 2),
			      sat.file = "results/04_Genome_wide_comparison_to_nonB_motifs/hsat_2_3_annotations/Pabelii.bed",
			      sat.title = "hsat2,3")	

ggsave("plots/Figure_S2_Genome_wide_plots/Pabelii.svg", P.Pabelii, width = 8, height = 2.5, units = "in", scale = 2.5)
	

#S. syndactylus

df = df.windows %>% filter(Species == "Ssyndactylus")

P.Ssyndactylus = Plot_reads_v_nonb(species = "Ssyndactylus",
			      cen = "resources/centromeres/Ssyndactylus.bed",
			      v.levels = c("Ssyndactylus_JAMBI_PDALseq", "Ssyndactylus_Karenina_PDALseq", "Ssyndactylus_JAMBI_TruSeq"),
			      v.labels = c("Ssyndactylus_JAMBI_PDALseq", "Ssyndactylus_Karenina_PDALseq", "Ssyndactylus_JAMBI_TruSeq"),
			      v.heights = c(0.25, 0.25, 1.2, 3),
			      sat.file = "results/04_Genome_wide_comparison_to_nonB_motifs/Ssyndactylus_subTerm_monhorSF4.bed",
			      sat.title = "subterminal alpha sat")	

ggsave("plots/Figure_S2_Genome_wide_plots/Ssyndactylus.svg", P.Ssyndactylus, width = 8, height = 3, units = "in", scale = 2.5)
	


