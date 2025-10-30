
library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)
library(tidyverse)
library(viridis)
library(cowplot)

####HSat3 synteny####

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


x = "hsa9"

P.A = ggplot() +
		ggtitle(x) +
		facet_wrap(~Species, nrow = 1) +
		geom_rect(data = df_size %>% filter(hsa == x), mapping = aes(ymax = 0, ymin = -end, xmax = 1.4, xmin = 1.1), color = "black", fill = "white") +
		geom_text(data = df_size %>% filter(hsa == x), mapping = aes(y = -0.75*end, x = 1.25, label = chr), color = "black", angle = 90, size = 2) +
		geom_rect(data = df_cen %>% filter(hsa == x), mapping = aes(ymax = -start, ymin = -end, xmin = 1.0, xmax = 1.5), fill = "dimgrey", color = "dimgrey") +
		geom_rect(data = df_Hsat23 %>% filter(hsa == x), mapping = aes(ymax = -start, ymin = -end, xmax = 0.95, xmin = 0.55), fill = "black") +
		geom_rect(data = df_DR %>% filter(hsa == x), mapping = aes(ymax = -start, ymin = -end, xmin = 0.05, xmax = 0.45, fill = density)) +
		geom_rect(data = df_PDAL %>% filter(hsa == x), mapping = aes(ymax = -start, ymin = -end, xmax = 0, xmin = -Norm_reads, fill = Norm_reads)) +
		theme_void() +
		theme(legend.position = "none") +
		scale_fill_gradient(low = "white", high = "Red") +
		theme(plot.title = element_text(size = 8),
		      strip.text = element_text(size = 6))
P.A 

####Hsat monomer count####

T2T = c("Hsapien_GM24385", "Ptroglodytes_AG18354", "Ppaniscus_KB8711", "Ggorilla_KB3781", "Ppygmaeus_AG05252", "Pabelii_AG06213", "Ssyndactylus_JAMBI")
T2T.labels = c("H. sapiens\nGM24385", "P. troglodytes\nAG18354", "P. paniscus\nKB8711", "G. gorilla\nKB3781", "P. pygmaeus\nAG05252", "P. abelii\nAG06213", "S. syndactylus\n JAMBI")

df.sat.count = read.csv("results/06_reference_free_repeat_analysis/Short_tandem_repeats_by_technical_replicate.csv") %>% filter(Repeat == "AATGG") 
df.sat.count = df.sat.count %>% filter(Bio_sample %in% T2T)

df.sat.count$Bio_sample = factor(df.sat.count$Bio_sample, levels = T2T, labels = T2T.labels) 

unique(paste(df.sat.count$Bio_sample, df.sat.count$Technique, sep = " "))

df.sat.count$Experiment = factor(paste(df.sat.count$Bio_sample, df.sat.count$Technique, sep = " "), levels = c("H. sapiens\nGM24385 PDAL-seq", "H. sapiens\nGM24385 TruSeq", "P. troglodytes\nAG18354 PDAL-seq", "P. troglodytes\nAG18354 TruSeq", "P. paniscus\nKB8711 PDAL-seq", "P. paniscus\nKB8711 TruSeq", "G. gorilla\nKB3781 PDAL-seq", "G. gorilla\nKB3781 TruSeq", "P. pygmaeus\nAG05252 PDAL-seq", "P. pygmaeus\nAG05252 TruSeq", "P. abelii\nAG06213 PDAL-seq", "P. abelii\nAG06213 TruSeq", "S. syndactylus\n JAMBI PDAL-seq", "S. syndactylus\n JAMBI TruSeq"), labels = c("H. sapiens PDAL-seq", "H. sapiens TruSeq", "P. troglodytes PDAL-seq", "P. troglodytes TruSeq", "P. paniscus PDAL-seq", "P. paniscus TruSeq", "G. gorilla PDAL-seq", "G. gorilla TruSeq", "P. pygmaeus PDAL-seq", "P. pygmaeus TruSeq", "P. abelii PDAL-seq", "P. abelii TruSeq", "S. syndactylus PDAL-seq", "S. syndactylus TruSeq"))

P.sat.count = ggplot(df.sat.count , aes(y = kbp_per_Mbp, x = Experiment, fill = Technique, group = Technique)) +
    geom_bar(stat="summary", fun = mean, color = "black") +
    stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
    geom_jitter(size = 0.5) +
    theme_classic() +
    theme() +
    ylab("Sequencing rate (kbp/Mbp)") + 
    ggtitle("Reference free HSat3 monomers") +
    scale_fill_manual(values = c("white", "grey")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    	  axis.title.x = element_blank(),
    	  axis.text.y = element_text(size = 6),
    	  axis.title.y = element_text(size = 8),
    	  legend.position= c(0.75,0.75),
    	  legend.title = element_text(size = 6),
    	  legend.text = element_text(size = 6),
    	  legend.key.size = unit(0.1, "in"),
    	  plot.title =  element_text(size = 8)) 
P.sat.count

###CD plots####


df = read.csv("results/18_biophysical_charaterization_of_HSat3_monomers/XX_Data_for_manuscript/CD_data/js3020_tidy_formated_CD_data.csv")

#Compare to published data

df.Chen = read.csv("results/18_biophysical_charaterization_of_HSat3_monomers/XX_Data_for_manuscript/CD_data/Chen_et_al_NAR_2025.csv")
colnames(df.Chen) = c("Wavelength_nm", "Molar_elipicity", "Buffer")
df.Chen$Study = "Yatsunyk et. al. NAR 2025 (ATGGA)4"

df.A = df %>% filter(Sample == "(AATGG)4") %>% select(Wavelength_nm, Molar_elipicity, Buffer)
df.A$Study = "This study (AATGG)4"
df.A = bind_rows(df.A, df.Chen)

PC = ggplot(df.A, aes(x = Wavelength_nm, y = Molar_elipicity, linetype = Buffer, color = Study)) +
  geom_line() +
  theme_classic() +
  ggtitle("CD of (AATGG)4 and (ATGGA)4\nin different buffer conditions") +
  ylab("Molar ellipticity (1/M*cm)") +
  xlab("Wavelength (nm)") +
  scale_x_continuous(breaks = c(10*22:33), labels = c("220", "", "240", "", "260", "", "280", "", "300", "", "320", ""))+
  scale_color_manual(values = inferno(3)) +
  scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    	  axis.title.x = element_text(size = 8),
    	  axis.text.y = element_text(size = 6),
    	  axis.title.y = element_text(size = 8),
    	  legend.position = "bottom",
    	  legend.title = element_text(size = 6),
    	  legend.text = element_text(size = 6),
    	  legend.key.size = unit(0.2, "in"),
    	  plot.title =  element_text(size = 8)) 
PC 

df = df %>% filter(Sample != "Buffer")

df$Buffer = factor(df$Buffer, levels = c("140 mM LiCl", "100 mM KCl 140 mM LiCl"))

PD = ggplot(df, aes(x = Wavelength_nm, y = CD_background_subtracted, linetype = Buffer, color = Sample)) +
 geom_line() +
  theme_classic() +
  scale_color_manual(values = magma(6)) +
  ggtitle("Effect of (AATGG)n repeat length\n and KCl on structure") +
  ylab("CD (mDeg)") +
  xlab("Wavelength (nm)") +
  scale_x_continuous(breaks = c(10*22:33), labels = c("220", "", "240", "", "260", "", "280", "", "300", "", "320", "")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    	  axis.title.x = element_text(size = 8),
    	  axis.text.y = element_text(size = 6),
    	  axis.title.y = element_text(size = 8),
    	  legend.position = "bottom",
    	  legend.title = element_text(size = 6),
    	  legend.text = element_text(size = 6),
    	  legend.key.size = unit(0.2, "in"),
    	  plot.title =  element_text(size = 8)) 
PD 

####Calculate Concentration for each sample####

####PCA####


Experiments = unique(df$Experiment)
df.PCA = df %>% filter(Sample != "Buffer")
l.df.PCA = df.PCA %>% group_by(Experiment) %>% group_split()
M.PCA = matrix(ncol = nrow(l.df.PCA[[1]]), nrow = length(Experiments))
colnames(M.PCA) = l.df.PCA[[1]]$Wavelength_nm
df$Experiment = paste(df$Sample, df$Buffer, sep = "_")
Experiments = c()
for(i in 1:length(l.df.PCA)){
  M.PCA[i,] = l.df.PCA[[i]]$CD_background_subtracted
  Experiments[i] = l.df.PCA[[i]]$Experiment[1]
}

rownames(M.PCA) = Experiments

results <- prcomp(M.PCA)
results$x
results$sdev^2 / sum(results$sdev^2)

df.PCA.Result = data.frame(results$x)
l.Experiment = strsplit(rownames(df.PCA.Result), split = "_")
Sample = c()
Buffer = c()
for (i in 1:length(l.Experiment)){
  Sample[i] = l.Experiment[[i]][1]
  Buffer[i] = l.Experiment[[i]][2]
}
df.PCA.Result$Sample = Sample
df.PCA.Result$Buffer = factor(Buffer, levels = c("140 mM LiCl", "100 mM KCl 140 mM LiCl"))

PE = ggplot(df.PCA.Result, aes(x = PC1, y = PC2, color = Sample, shape = Buffer)) +
  geom_point(size = 1.0) +
  scale_color_manual(values =magma(6)) +
  theme_classic() +
  xlab("PC1 Repeat Length") +
  ylab("PC2 Odd/Even") +
  ggtitle("Principle components analysis") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    	  axis.title.x = element_text(size = 8),
    	  axis.text.y = element_text(size = 6),
    	  axis.title.y = element_text(size = 8),
    	  legend.position = "bottom",
    	  legend.title = element_text(size = 6),
    	  legend.text = element_text(size = 6),
    	  legend.key.size = unit(0.2, "in"),
    	  plot.title =  element_text(size = 8)) 
PE  

#####SEC MALS####

df.SEC = read.csv("results/18_biophysical_charaterization_of_HSat3_monomers/js3022_SEC_UV_chromatograms_tidy_formatted.csv")
P.chrom = ggplot() +
  geom_rect(data = df.SEC, mapping = aes(ymin = Min, ymax = Max, xmin = -0.5, xmax = -0.1, fill = UV))  +
  geom_point(data = df.SEC, mapping = aes(y = Min, x  = UV), color = "black", size = 0.25)  +
  theme_classic() + 
  facet_wrap(~Sample, nrow = 1) +
  scale_fill_gradient(low = "white", high = "black") +
  scale_y_reverse(limits = c(45, 20)) +
  ggtitle("Size Exclusion Chromatography (SEC) - UV absorbance") +
  xlab("Absorbance (260 nm)") +
  ylab("Retention time (min)")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    	  axis.title.x = element_text(size = 8),
    	  axis.text.y = element_text(size = 6),
    	  axis.title.y = element_text(size = 8),
    	  legend.position = "bottom",
    	  legend.title = element_text(size = 6),
    	  legend.text = element_text(size = 6),
    	  legend.key.size = unit(0.1, "in"),
    	  plot.title =  element_text(size = 8),
    	  strip.text = element_text(size = 7)) 
P.chrom


####Compile plots####

PAB = plot_grid(P.A, P.sat.count, nrow = 1, labels = c("A", "B"), rel_widths = c(1.5,1), label_size = 12)
PCDE = plot_grid(PC, PD, PE, nrow = 1, labels = c("C", "D", "E"), label_size = 12, align = "v")
PFG = plot_grid(P.chrom, P.chrom, nrow = 1, labels = c("F", "G"), rel_widths = c(1.4,1), label_size = 12)
P = plot_grid(PAB, PCDE, PFG, rel_heights = c(1,1,1.2), ncol = 1)

ggsave("new_plots/Figure_6.svg", P, width = 6.5, height = 8, units = "in")



