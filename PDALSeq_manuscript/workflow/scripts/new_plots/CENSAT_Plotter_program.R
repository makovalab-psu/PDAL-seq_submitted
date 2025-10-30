
library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)
library(tidyverse)
library(viridis)
library(cowplot)

Species = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Pabelii", "Ppygmaeus")
Species_labels = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. abelii", "P. pygmaeus")


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

Exp.labels   = c("H. sapiens GM24385 PDAL-seq",
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
		 
Censat.levels = c("active_HOR", "inactive_HOR", "diverged_monomeric_HOR", "bSat", "HSat1", "HSat2", "HSat3", "Other_satellite", "rDNA", "Transition_region")
 Censat.labels = c("active HOR", "inactive HOR", "diverged/monomeric HOR", "bSat", "HSat1", "HSat2", "HSat3", "Other satellite", "rDNA", "Transition region")
Censat.colors = c("#891641", "#CF6732", "#F8DB5B", "#DCAFD0", "#18988B", "#323467", "#77A1BA", "#768661", "#000000", "#959595")
####Censat plot####

df.Censat = read.csv("results/17_PDALseq_reads_in_CenSat/Human_Cen_Sat_coordinates.csv")
df.Censat$Altemose = factor(df.Censat$Altemose, levels = Censat.levels, labels = Censat.labels)
df.chr.sizes = read.delim("resources/genomes/Hsapien.chrom.sizes", header = F)
colnames(df.chr.sizes) = c("chr", "end")
df.PDALseq = read.delim("data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/10000_nucleotides/Hsapien_GM24385_PDALseq.bg", header = F)
colnames(df.PDALseq) = c("chr", "start", "end", "reads")
Q99 = quantile(df.PDALseq$reads, 0.999)
Norm_reads = df.PDALseq$reads/Q99
Norm_reads[which(Norm_reads >1)] = 1
df.PDALseq$Norm_reads = Norm_reads

plot_PDALseq_cen = function(CHR = "chr9"){
	df.x = df.Censat  %>% filter(chr == CHR)
	Max_x = max(df.x %>% select(end)) + 1000000
	Min_x = min(df.x %>% select(end)) - 1000000
	x_censat = unique(df.x$Altemose)
	x_colors = Censat.colors[which(Censat.labels %in% x_censat)]
	P = ggplot() +
		ggtitle(CHR) + 
		geom_rect(mapping = aes(xmin = Min_x, xmax = Max_x, ymin = -0.25, ymax = 0), fill = "white", color = "black") +
		geom_rect(data = df.x, mapping = aes(xmin = start, xmax = end, ymin = -0.25, ymax = 0, fill = Altemose)) +
		geom_rect(data = df.PDALseq  %>% filter(chr == CHR) %>% filter(start >= Min_x) %>% filter(end <= Max_x), mapping = aes(xmin = start, xmax = end, ymin = 0, ymax = Norm_reads), fill = "red") +
		theme_void() +
		scale_fill_manual(values = x_colors) +
		theme(plot.title = element_text(size=8, color="black"),
		      legend.position = "bottom",
		      legend.text = element_text(size=6, color="black"),
		      legend.title = element_text(size=6, color="black")) 
		plot_path = paste("new_plots/Cen_sat_vrs_PDALseq/", CHR, ".png", sep = "")
		ggsave(plot_path, P, width = 6.5, height = 3, units = "in")
		return(P)
}

lapply(unique(df.Censat$chr), plot_PDALseq_cen)


plot_PDALseq_cen("chr22")





