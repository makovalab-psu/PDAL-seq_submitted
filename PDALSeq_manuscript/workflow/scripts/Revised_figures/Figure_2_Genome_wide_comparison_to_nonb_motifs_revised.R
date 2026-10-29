
library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)
v.margins = c(0.1,0.1,0.1,0.1)

theme_set(theme_classic(base_size = 6) +
            theme(text = element_text(color = "black", size = 6),
                  axis.text = element_text(color = "black", size = 6),
                  axis.title = element_text(color = "black", size = 6),
                  strip.text = element_text(color = "black", size = 6),
                  legend.text = element_text(color = "black", size = 6),
                  legend.title = element_text(color = "black", size = 6),
                  plot.title = element_text(color = "black", size = 6)))


####99% Winsorization####

df.windows = readr::read_csv("results/21_add_K562_to_analysis/nonB_coverage_versus_read.csv") %>% filter(nonB == "all")

unique(df.windows$chr)
df.windows$Experiment = paste(df.windows$Species, df.windows$Bio_sample, df.windows$Technique, sep = "_")
unique(df.windows$Experiment)
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
v.levels = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq", "Hsapien_K562_PDALseq")
v.labels = c("GM24385-HG002",  "GM24143", "HEK293T", "CCL86-Raji", "K562")


df.reads = df.windows %>% filter(Species == species)
df.reads$Experiment = factor(df.reads$Experiment, levels = v.levels, labels = v.labels) 

#Plot reads v coordinate
P.reads = ggplot(df.reads %>% filter(Experiment == "GM24385-HG002"), aes(x = interval, y = reads, fill = chr)) +
       facet_wrap(~Experiment, ncol = 1) +
	geom_bar(stat = "identity") +
	theme_classic() +
	scale_fill_manual(values = viridis(length(unique(df.reads$chr)), option = "turbo")) +
	scale_x_continuous(limits = c(min(df.reads$interval), max(df.reads$interval)), expand = c(0, 0))  +
  	theme(legend.position = "none",
          axis.title.y = element_text(size = 8, color = "black"),
          strip.text = element_text(size = 8, color = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank()) +
          ylab("PDAL-Seq reads")
	
#Plot non-B v coordinate
	
df.nonb = read.csv("results/21_add_K562_to_analysis/nonB_coverage_versus_read.csv") 
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

NONB = c("all", "GQ", "SLS", "DR", "STR", "Z", "TRI", "MR", "CRU", "IR", "APR")

NONB_labels = c("all", "G4", "SLS", "DR", "STR", "Z", "TRI", "MR", "CRU", "IR", "APR")

df.nonb = df.nonb %>% filter(nonB %in% NONB)
df.nonb$nonB = factor(df.nonb$nonB, levels = rev(NONB), labels = rev(NONB_labels))

P.motif = ggplot(df.nonb , aes(x = interval, y = nonB, fill = nonb.count)) +
        geom_tile() +
	theme_classic() +
	scale_fill_gradient(low = "white", high = "red") +
	scale_x_continuous(limits = c(min(df.nonb$interval), max(df.nonb$interval)), expand = c(0, 0)) +
    	theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.text = element_text(size = 6, color = "black"),
          legend.title = element_text(size = 6, color = "black"),
          legend.key.size = unit(0.5, "lines")) +
          ylab("non-B motifs")
	
	#Label in centromere coordinates
	df.cen = read.delim(cen, header = F)
	colnames(df.cen) = c("chr", "start", "end", "annotation")
	
l.df.nonb = df.nonb %>% group_by(chr) %>% group_split()
          
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
	
	
#Make an ideogram
P.ideogram = ggplot(df.centromere, aes(x = interval, y = 1, fill = chr)) +
        geom_tile() +
        geom_point(mapping = aes(x = interval, y = Centromere), size = 2, color = "grey") + 
	theme_classic() +
	scale_fill_manual(values = viridis(length(unique(df.centromere$chr)), option = "turbo")) +
	scale_x_continuous(limits = c(min(df.centromere$interval), max(df.centromere$interval)), expand = c(0, 0)) +
    	theme_void() +
    	theme(legend.position = "none",
    	      plot.title = element_text(size=8, color = "black")) 
    	


#Corrplot all, IR, TRI, G4#

Cor_nonB = c("all", "G4", "TRI", "CRU")

df.B = df.nonb %>% filter(nonB %in% Cor_nonB)
df.B$nonB = factor(df.B$nonB, levels = Cor_nonB)
unique(df.B$chr)
df.B$Sex_chr = factor(df.B$chr, levels = unique(df.B$chr), labels = c(rep("FALSE", 22), "chrX", "chrY", "FALSE"))


P.B = ggplot(df.B, aes(x = nonb.count, y = reads, color = Sex_chr)) +
	facet_wrap(~nonB) +
	theme_classic() +
	geom_point(alpha = 0.1, size = 0.5) +
	xlab("non-B DNA motif density") +
	ylab("PDAL-Seq read density") +
	scale_color_manual(values = viridis(4, option= "inferno"))  + 
  	geom_smooth(data = df.B %>% filter(Sex_chr == "FALSE"),
  	      method = "loess", 
              span = 10,
              se = FALSE,
              color = "red",
              size = 0.5) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size =6, color = "black"),
    	  axis.text.y = element_text(size = 6, color = "black"),
    	  axis.title.x = element_text(size = 8, color = "black"),
    	  axis.title.y = element_text(size = 8, color = "black"),
    	  strip.text = element_text(size = 8, color = "black"),
          legend.position = "bottom",
          legend.text = element_text(size = 8, color = "black"),
          legend.title = element_blank())

####Genome wide correlation####

df = read.csv("results/21_add_K562_to_analysis/nonB_correlation_coefficient.csv") %>% filter(windows == 1000000)

df$Species = factor(df$Species, levels = c("Hsapien","Ptroglodytes", "Ppaniscus",    "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus"), labels = c("H. sapiens","P. troglodytes", "P. paniscus",    "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus")) 
df$nonB = factor(df$nonB, levels = rev(c("APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS",  "GQ", "all")), labels = rev(c("APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS",  "G4", "all")))
df$Technique = factor(df$Technique, levels = c("PDALseq", "TruSeq"), labels = c("PDAL-Seq", "Control (TruSeq)"))

Species_palette = c("#9E9E9E", "#BE242F","#DC93B6", "#1C7D90", "#F2D983", "#FA9632", "#313757")

P.cor = ggplot(data = df, mapping = aes(x = nonB, y = Spearman)) +
	facet_wrap(~Technique, nrow = 1) +
	geom_bar(stat = "summary", fun = mean, fill = "grey") +
	stat_summary(fun.data = mean_sd, geom = "errorbar", color = "black") +
	geom_beeswarm(mapping = aes(color = Species), size = 0.4) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, size =6, color = "black"),
    	  axis.text.y = element_text(size = 6, color = "black"),
    	  axis.title.x = element_blank(),
    	  axis.title.y = element_text(size = 8, color = "black"),
    	  strip.text = element_text(size = 8, color = "black"),
          legend.position = "none",
          legend.text = element_text(size = 8, color = "black"),
          legend.title = element_blank()) +
	scale_color_manual(values = Species_palette) +
	geom_hline(yintercept = 0)      +
	ylab("Spearman's rank\ncorrelation")   

####Fraction in syntenic versus non-syntenic####

df = readr::read_csv("results/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/PDAL_Seq_signal_in_and_out_of_syntenic_regions.csv")

head(df)
df$fraction_of_total_reads = df$Bases_in_region/(df$Total_bases)
df$Normalized_reads = df$Bases_in_region/(df$Total_bases*df$Region_nt)
unique(df$Region)
df$Region = factor(df$Region, levels = c("promoter", "utr5", "exon", "intron", "utr3", "censat", "nongenic_repeats", "nongenic_nonrepetitive"))

df$Method = NA
df$Method[which(grepl("PDALseq", df$Dataset))] = "PDAL-Seq"
df$Method[which(grepl("TruSeq", df$Dataset))] = "TruSeq"

Experiments = unique(df$Dataset)

l_df = {}
for (i in 1:length(Experiments)){
  l_df[[i]] = data.frame(df %>% filter(Dataset == Experiments[i]) %>% group_by(Syntenic) %>% summarise(sum = sum(fraction_of_total_reads)))
  l_df[[i]]$Experiment = Experiments[i]
  l_df[[i]]$Species = strsplit(Experiments[i], split = "_")[[1]][1] 
}
df_sum = bind_rows(l_df)
df_sum$Method = NA
df_sum$Method[which(grepl("PDALseq", df_sum$Experiment))] = "PDAL-Seq"
df_sum$Method[which(grepl("TruSeq", df_sum$Experiment))] = "TruSeq"
df_sum$Syntenic = factor(df_sum$Syntenic, levels = c("Syntenic", "Non-syntenic"))

df_sum$Species = factor(df_sum$Species, levels = c("Hsapien","Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus"), labels = c("H. sapiens","P. troglodytes", "P. paniscus", "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus"))
Species_palette = c("#9E9E9E", "#BE242F","#DC93B6", "#1C7D90", "#F2D983", "#FA9632", "#313757")

PA = ggplot(df_sum, aes(x = Syntenic, y = sum)) +
  geom_bar(stat = "summary", fun = "mean") +
  stat_summary(fun.data = mean_sd, geom = "errorbar", color = "black") +
  geom_jitter(mapping = aes(color = Species), size = 0.4) +
  facet_wrap(~Method, nrow = 1) + 
  scale_color_manual(values = Species_palette) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0,0.2,0.4,0.6, 0.8, 1.0)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  ylab("Fraction of reads")


####Compile plot####
v.heights = c(0.1, 1.0, 1.4)

P.A = plot_grid(P.ideogram, P.reads,P.motif, align = "v", ncol = 1, rel_heights = v.heights, labels = c("A"), label_size = 8)
P.BC = plot_grid(PA, P.cor, P.B, nrow = 1, labels = c("B", "C", "D"), label_size = 8, rel_widths = c(0.3,0.8,1), align = "h")
P = plot_grid(P.A, P.BC, ncol = 1, rel_heights = c(1.2,1), label_size = 8)          

ggsave("revised_plots/Figure_2.svg", P, width = 6.5, height = 4.5, units = "in")


