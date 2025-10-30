
library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)
library(tidyverse)
library(viridis)
library(cowplot)

Species = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Pabelii", "Ppygmaeus", "Ssyndactylus")
Species_labels = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. abelii", "P. pygmaeus", "S. syndactylus")
Species_palette = c("#9E9E9E", "#BE242F","#DC93B6", "#1C7D90", "#F2D983", "#FA9632", "#313757")


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
Censat.labels = c("active HOR", "inactive HOR", "div/mon HOR", "bSat", "HSat1", "HSat2", "HSat3", "Other satellite", "rDNA", "Transition region")
Censat_y =  c("active HOR", "inactive HOR", "div/mon HOR", "bSat", "HSat1", "HSat2", "HSat3", "Other satellite", "Transition region", "rDNA")
Censat.colors = c("#891641", "#CF6732", "#F8DB5B", "#DCAFD0", "#18988B", "#323467", "#77A1BA", "#768661", "#000000", "#959595")

nonB.levels = c("GQ", "SLS", "DR", "STR", "Z", "TRI", "MR", "CRU", "IR", "APR")
nonB.labels = c("G4", "SLS", "DR", "STR", "Z", "TRI", "MR", "CRU", "IR", "APR")

####Censat plot####

####Sat array enrichment####

df.E = read.csv("results/17_PDALseq_reads_in_CenSat/08_Altemose_paper_censat_enrichment_classes.csv")
df.E$Experiment = factor(df.E$Experiment, levels = Exp.levels, labels = Exp.labels)
df.E$Censat = factor(df.E$Censat, levels = rev(Censat.levels), labels = rev(Censat.labels))
df.E = df.E %>% filter(!is.na(Enrichment))
df.E$log2Enrichment[which(df.E$log2Enrichment > log2(3))] = log2(3)
df.E$log2Enrichment[which(df.E$log2Enrichment < -log2(3))] = -log2(3)
Sig = rep("ns", nrow(df.E))
Sig[which(df.E$Wilcox.p <= 0.05)] = "*"
Sig[which(df.E$Wilcox.p <= 0.01)] = "**"
df.E$Sig.labels = paste(round(df.E$Enrichment, digits = 1), Sig, sep = "\n")

Max.E = max(abs(df.E %>% filter(Censat != "rDNA") %>% select(log2Enrichment)))

df.PDAL =  df.E %>% filter(Experiment %in% c("H. sapiens GM24385 PDAL-seq", "H. sapiens GM24143 PDAL-seq", "P. troglodytes AG18354 PDAL-seq", "P. troglodytes AG18358 PDAL-seq", "P. paniscus KB8711 PDAL-seq", "P. paniscus KB14571 PDAL-seq", "G. gorilla KB3781 PDAL-seq", "G. gorilla F7549 PDAL-seq", "P. pygmaeus AG05252 PDAL-seq", "P. abelii AG06213 PDAL-seq", "S. syndactylus JAMBI PDAL-seq", "S. syndactylus Karenina PDAL-seq")) #%>% filter(Censat != "rDNA")

df.Tru =  df.E %>% filter(Experiment %in% c("H. sapiens GM24385 TruSeq", "P. troglodytes AG18354 TruSeq", "P. paniscus KB8711 TruSeq", "G. gorilla KB3781 TruSeq", "P. pygmaeus AG05252 TruSeq", "P. abelii AG06213 TruSeq", "S. syndactylus JAMBI TruSeq")) #%>% filter(Censat != "rDNA")

df.int = read.delim("results/17_PDALseq_reads_in_CenSat/Altemose_paper_Censat_nonB_intersect.txt") #%>% filter(censat != "rDNA")
df.int = df.int[-which(is.na(df.int$total)),]
df.int$intersect[which(is.na(df.int$intersect))] = 0
df.int$Percent =  100*df.int$intersect/df.int$total
df.int$censat = factor(df.int$censat, levels = Censat.levels, labels = Censat.labels)
df.int$nonb = factor(df.int$nonb, levels = nonB.levels, labels = nonB.labels) 
df.int = df.int %>% filter(nonb %in% c("G4", "DR", "Z", "TRI", "CRU", "APR"))

P.Tru  = ggplot(df.Tru, aes(x = Experiment, y = Censat, fill = log2Enrichment, label = Sig.labels )) +
	geom_tile() +
	geom_text(size = 2) +
	theme_classic() +
	scale_y_discrete(limits = rev(Censat_y)) +
	scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-log2(3), log2(3))) +
              theme(axis.text.x = element_text(size=6, color="black", angle=45, hjust = 1),
	      axis.title.x = element_blank(),
	      axis.text.y = element_text(size=6, color="black"),
	      axis.title.y = element_blank(),
	      plot.title = element_text(size=8, color="black"),
	      legend.position = "top",
	      legend.key.size = unit(0.2, "in"),
    	      legend.title = element_text(size=8, color="black"),
	      legend.text = element_text(size=8, color="black"))

P.PDAL = ggplot(df.PDAL, aes(x = Experiment, y = Censat, fill = log2Enrichment, label = Sig.labels )) +
	geom_tile() +
	geom_text(size = 2) +
	theme_classic() +
	scale_y_discrete(limits = rev(Censat_y)) +
	scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-log2(3), log2(3))) + 
              theme(axis.text.x = element_text(size=6, color="black", angle=45, hjust = 1),
	      axis.title.x = element_blank(),
	      axis.text.y = element_text(size=6, color="black"),
	      axis.title.y = element_blank(),
	      plot.title = element_text(size=8, color="black"),
	      legend.key.size = unit(0.1, "in"),
    	      legend.position = "top",
	      legend.title = element_text(size=8, color="black"),
	      legend.text = element_text(size=6, color="black"),
	      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) 

df.SI = bind_rows(df.PDAL, df.Tru)


P.SI = ggplot(df.SI, aes(x = Biosample, y = Censat, fill = log2Enrichment, label = Sig.labels )) +
	facet_wrap(~Technique, ncol = 1) +
	geom_tile() +
	geom_text(size = 2) +
	theme_classic() +
	scale_x_discrete(limits = unique(df.SI$Biosample)) +
	scale_y_discrete(limits = rev(Censat_y)) +
	scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-log2(3), log2(3))) + 
              theme(axis.text.x = element_text(size=6, color="black", angle=45, hjust = 1),
	      axis.title.x = element_blank(),
	      axis.text.y = element_text(size=6, color="black"),
	      axis.title.y = element_blank(),
	      plot.title = element_text(size=8, color="black"),
	      legend.key.size = unit(0.1, "in"),
    	      legend.position = "top",
	      legend.title = element_text(size=8, color="black"),
	      legend.text = element_text(size=6, color="black"),
	      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) 


Paired_samples = c("GM24385",
		 "AG18354",
		 "KB8711",
		 "KB3781",
		 "AG05252",
		 "AG06213",
		 "JAMBI")


l.df = {}
p.Wilcoxon_signed_rank_test = c()
for (i in 1:length(Censat_y)){
	Censat.i = Censat_y[i]
	df.PDAL.i = df.PDAL %>% filter(Biosample %in% Paired_samples) %>% filter(Censat == Censat.i) %>% arrange(Biosample)
	df.Tru.i = df.Tru %>% filter(Biosample %in% Paired_samples) %>% filter(Censat == Censat.i) %>% arrange(Biosample)
	PDALseq_enrichment = df.PDAL.i$Enrichment
	Test.wilcox = NA
	TruSeq_enrichment = df.Tru.i$Enrichment
	Test.wilcox = wilcox.test(PDALseq_enrichment, TruSeq_enrichment, paired=TRUE)
	print(Test.wilcox)
	p.Wilcoxon_signed_rank_test[i] = Test.wilcox$p.value
	df.i = df.PDAL.i[,c(1:4)]
	df.i$PDALseq_enrichment = PDALseq_enrichment
	df.i$TruSeq_enrichment = TruSeq_enrichment
	df.i$Increase_in_PDALseq = PDALseq_enrichment/TruSeq_enrichment
	df.i$Censat = Censat.i 
	l.df[[i]] = df.i
}

df.compare = bind_rows(l.df)
df.compare$log2_Increase_in_PDALseq = log2(df.compare$Increase_in_PDALseq)
df.compare$log2_Increase_in_PDALseq[which(df.compare$log2_Increase_in_PDALseq > log2(3))] = log2(3)
df.compare$log2_Increase_in_PDALseq[which(df.compare$log2_Increase_in_PDALseq < -log2(3))] = -log2(3)
df.compare$En_labels = as.character(round(df.compare$Increase_in_PDALseq, digits = 2))

P.compare = ggplot(df.compare, aes(x = Experiment, y = Censat, fill = log2_Increase_in_PDALseq, label = En_labels )) +
	geom_tile() +
	geom_text(size = 2) +
	theme_classic() +
	scale_y_discrete(limits = rev(Censat_y)) +
	scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-log2(3), log2(3))) + 
              theme(axis.text.x = element_text(size=6, color="black", angle=45, hjust = 1),
	      axis.title.x = element_blank(),
	      axis.text.y = element_text(size=6, color="black"),
	      axis.title.y = element_blank(),
	      plot.title = element_text(size=8, color="black"),
	      legend.key.size = unit(0.1, "in"),
    	      legend.position = "top",
	      legend.title = element_text(size=8, color="black"),
	      legend.text = element_text(size=6, color="black"),
	      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) 

df.sig = data.frame("Censat" = Censat_y, "p" = p.Wilcoxon_signed_rank_test, "Sig.label" = paste("p =", round(p.Wilcoxon_signed_rank_test, digits = 3)))

P.sig = ggplot(df.sig, aes(x = 1, y = Censat, label =  Sig.label)) +
	geom_text(size = 2) +
	theme_void() +
	scale_y_discrete(limits = rev(Censat_y)) +
	theme(plot.margin = margin(t = 0.5, r = 0, b = 0.75, l = 0, unit = "in"))



df.int$species = factor(df.int$species, level = Species)	      
P.int = ggplot(df.int, aes(x = Percent, y = censat)) +
    geom_bar(stat="summary", fun = mean, color = "black") +
    stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
    geom_beeswarm(mapping = aes(color = species),size = 0.5) +
    theme_classic() +
    facet_wrap(~nonb, nrow = 1, scale = "free_x") +
    scale_color_manual(values = Species_palette) +
    scale_y_discrete(limits = rev(Censat_y)) +
              theme(axis.text.x = element_text(size=6, color="black", angle=90),
	      axis.title.x = element_text(size=8, color="black"),
	      axis.text.y = element_blank(),
	      axis.title.y = element_blank(),
	      plot.title = element_text(size=8, color="black"),
	      legend.position = "none",
	      legend.title = element_text(size=8, color="black"),
	      legend.text = element_text(size=8, color="black"),
	      strip.text = element_text(size=8, color="black"),
	      plot.margin = margin(t = 0.325, r = 0, b = 0.85, l = 0, unit = "in")) +
    xlab("Percent non-B DNA motif in satallite (%)")

P.A = plot_grid(P.compare, P.sig, P.int, align = "v", rel_widths = c(6,0.5,6), nrow = 1) 

    

####Hsat monomer count####

T2T = c("Hsapien_GM24385", "Ptroglodytes_AG18354", "Ppaniscus_KB8711", "Ggorilla_KB3781", "Ppygmaeus_AG05252", "Pabelii_AG06213", "Ssyndactylus_JAMBI")
T2T.labels = c("H. sapiens\nGM24385", "P. troglodytes\nAG18354", "P. paniscus\nKB8711", "G. gorilla\nKB3781", "P. pygmaeus\nAG05252", "P. abelii\nAG06213", "S. syndactylus\n JAMBI")

df.sat.count = read.csv("results/06_reference_free_repeat_analysis/Short_tandem_repeats_by_technical_replicate.csv") %>% filter(Repeat == "AATGG") 
df.sat.count = df.sat.count %>% filter(Bio_sample %in% T2T)

df.sat.count$Bio_sample = factor(df.sat.count$Bio_sample, levels = T2T, labels = T2T.labels) 

unique(paste(df.sat.count$Bio_sample, df.sat.count$Technique, sep = " "))

df.sat.count$Experiment = factor(paste(df.sat.count$Bio_sample, df.sat.count$Technique, sep = " "), levels = c("H. sapiens\nGM24385 PDAL-seq", "H. sapiens\nGM24385 TruSeq", "P. troglodytes\nAG18354 PDAL-seq", "P. troglodytes\nAG18354 TruSeq", "P. paniscus\nKB8711 PDAL-seq", "P. paniscus\nKB8711 TruSeq", "G. gorilla\nKB3781 PDAL-seq", "G. gorilla\nKB3781 TruSeq", "P. pygmaeus\nAG05252 PDAL-seq", "P. pygmaeus\nAG05252 TruSeq", "P. abelii\nAG06213 PDAL-seq", "P. abelii\nAG06213 TruSeq", "S. syndactylus\n JAMBI PDAL-seq", "S. syndactylus\n JAMBI TruSeq"), labels = c("H. sapiens PDAL-seq", "H. sapiens TruSeq", "P. troglodytes PDAL-seq", "P. troglodytes TruSeq", "P. paniscus PDAL-seq", "P. paniscus TruSeq", "G. gorilla PDAL-seq", "G. gorilla TruSeq", "P. pygmaeus PDAL-seq", "P. pygmaeus TruSeq", "P. abelii PDAL-seq", "P. abelii TruSeq", "S. syndactylus PDAL-seq", "S. syndactylus TruSeq"))

df.sat.count$Technique = factor(df.sat.count$Technique, levels = c("PDAL-seq", "TruSeq"), labels = c("PDAL-Seq", "Control (TruSeq)"))

P.sat.count = ggplot(df.sat.count, aes(y = kbp_per_Mbp, x = Experiment, fill = Technique, group = Technique)) +
    geom_bar(stat="summary", fun = mean, color = "black") +
    stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
    geom_jitter(size = 0.5) +
    theme_classic() +
    theme() +
    ylab("Sequencing rate\n(kbp/Mbp)") + 
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
    	  legend.position = "none",
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
    	  legend.position = "none",
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
    	  legend.position = "none",
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
    	  legend.position = "right",
    	  legend.title = element_text(size = 6),
    	  legend.text = element_text(size = 6),
    	  legend.key.size = unit(0.1, "in"),
    	  plot.title =  element_text(size = 8),
    	  strip.text = element_text(size = 7)) 
P.chrom


####Compile plots####
PA = plot_grid(P.compare, P.sig, P.int, align = "v", rel_widths = c(6,1,6), labels = c("A"), label_size = 12, nrow = 1)

PBCD = plot_grid(P.sat.count, PD, PE, nrow = 1, labels = c("B", "C", "D"), rel_widths = c(1,1,1), label_size = 12)

PEF = plot_grid(P.chrom, PE, nrow = 1, labels = c("E", "F"), rel_widths = c(2,1),  label_size = 12)
P = plot_grid(PA, PBCD, PEF, rel_heights = c(1.4,1,1), ncol = 1)

ggsave("new_plots/Figure_6.svg", P, width = 6.5, height = 8, units = "in")



