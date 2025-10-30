library(tidyverse)
library(viridis)
library(cowplot)

Species = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Pabelii", "Ppygmaeus")
Species_labels = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. abelii", "P. pygmaeus")

Experiment = c("Hsapien_GM24385_PDALseq", "Ptroglodytes_AG18354_PDALseq", "Ppaniscus_KB8711_PDALseq", "Ggorilla_KB3781_PDALseq", "Pabelii_AG06213_PDALseq", "Ppygmaeus_AG05252_PDALseq")



l.df_size = {}
l.df_cen = {}
l.df_Hsat23 = {}
l.df_PDAL = {} 
l.df_DR = {} 

for (i in 1:length(Species)){
	x = Species[i]
	#Chromosome sizes
	size_path = paste("resources/genomes/", x, ".chrom.sizes", sep = "") 
	df_size = readr::read_delim(size_path, col_names = c("chr", "end"))
	df_size[-which(grepl("utig",  df_size$chr)),]
	CHR = df_size$chr
	df_size$chr = factor(df_size$chr, levels = CHR)
	df_size$Species = x
	
	HSA = c()
	l.chr = strsplit(CHR, split = "_") 
	if (length(l.chr[[1]])==1){HSA = gsub("chr", "hsa", CHR)}else{
		for (j in 1:length(l.chr)){
			HSA[j] = l.chr[[j]][3]
		}
	}
	df_size$hsa = factor(df_size$chr, levels = CHR, labels = HSA)
	
	#Centromeres
	cen_path = paste("resources/centromeres/", x, ".bed", sep = "") 
	df_cen = readr::read_delim(cen_path, col_names = c("chr", "start", "end", "class"))
	df_cen$chr = factor(df_cen$chr, levels = CHR)
	df_cen$Species = x
	df_cen$hsa = factor(df_cen$chr, levels = CHR, labels = HSA)
	
	#Hsat23
	Hsat23_path = paste("results/04_Genome_wide_comparison_to_nonB_motifs/hsat_2_3_annotations/", x, ".bed", sep = "") 
	df_Hsat23 = readr::read_delim(Hsat23_path, col_names = c("chr", "start", "end"))
	df_Hsat23$chr = factor(df_Hsat23$chr, levels = CHR)
	df_Hsat23$Species = x
	df_Hsat23$hsa = factor(df_Hsat23$chr, levels = CHR, labels = HSA)
	
	#DR
	DR_path = paste("data/04_Genome_wide_comparison_to_nonB_motifs/03_nonb_density/1000000_nucleotides/", x, "/DR.bg", sep = "") 
	df_DR = readr::read_delim(DR_path, col_names = c("chr", "start", "end", "count", "bp", "length", "density"))
	df_DR$chr = factor(df_DR$chr, levels = CHR)
	df_DR$Species = x
	df_DR$hsa = factor(df_DR$chr, levels = CHR, labels = HSA)
	
	#PDAL-seq
	PDAL_path = paste("data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/1000000_nucleotides/", Experiment[i], ".bg", sep = "") 
	df_PDAL = readr::read_delim(PDAL_path, col_names = c("chr", "start", "end", "reads"))
	df_PDAL$chr = factor(df_PDAL$chr, levels = CHR)
	df_PDAL$Species = x
	df_PDAL$hsa = factor(df_PDAL$chr, levels = CHR, labels = HSA)
	df_PDAL$Experiment = Experiment[i]
	Q99 = quantile(df_PDAL$reads, 0.99)
	Norm_reads = df_PDAL$reads/Q99
	Norm_reads[which(Norm_reads > 1)] = 1
	df_PDAL$Norm_reads = Norm_reads

	l.df_size[[i]] = df_size
	l.df_cen[[i]] = df_cen
	l.df_Hsat23[[i]] = df_Hsat23
	l.df_DR[[i]] = df_DR
	l.df_PDAL[[i]] = df_PDAL
}

df_size = bind_rows(l.df_size)
df_cen = bind_rows(l.df_cen)
df_Hsat23 = bind_rows(l.df_Hsat23)
df_DR = bind_rows(l.df_DR)
df_PDAL = bind_rows(l.df_PDAL)

df_size$Species = factor(df_size$Species, levels = Species, labels = Species_labels)
df_cen$Species = factor(df_cen$Species, levels = Species, labels = Species_labels)
df_Hsat23$Species = factor(df_Hsat23$Species, levels = Species, labels = Species_labels)
df_DR$Species = factor(df_DR$Species, levels = Species, labels = Species_labels)
df_PDAL$Species = factor(df_PDAL$Species, levels = Species, labels = Species_labels)

plot_hsa = function(x = "hsa9"){
	P = ggplot() +
		ggtitle(x) +
		geom_rect(data = df_size %>% filter(hsa == x), mapping = aes(xmin = 0, xmax = end, ymin = -1.1, ymax = -1.9), color = "black", fill = "white") +
		geom_text(data = df_size %>% filter(hsa == x), mapping = aes(x = 0.75*end, y = -1.5, label = chr), color = "black") +
		geom_rect(data = df_cen %>% filter(hsa == x), mapping = aes(xmin = start, xmax = end, ymin = -1.0, ymax = -2.00), fill = "dimgrey", color = "dimgrey") +
		geom_rect(data = df_Hsat23 %>% filter(hsa == x), mapping = aes(xmin = start, xmax = end, ymin = -0.95, ymax = -0.55), fill = "black", color = "black") +
		geom_rect(data = df_DR %>% filter(hsa == x), mapping = aes(xmin = start, xmax = end, ymin = -0.45, ymax = -0.05, fill = density)) +
		geom_rect(data = df_PDAL %>% filter(hsa == x), mapping = aes(xmin = start, xmax = end, ymin = 0, ymax = Norm_reads, fill = Norm_reads)) +
		facet_wrap(~Species, ncol = 1) +
		theme_void() +
		scale_x_continuous(limits = c(0, 253847350)) +
		theme(legend.position = "none") +
		scale_fill_gradient(low = "white", high = "Red") 
	plot_path = paste("plots/Figure_SX_hsat_synteny/", x, ".svg", sep = "")
	ggsave(plot_path, P, width = 5, height = 6)	
	return(P)
}

l.P = lapply(unique(df_size$hsa), plot_hsa)

