####Figure 1D Correlation####

library(tidyverse)
library(viridis)
library(cowplot)

args = commandArgs(trailingOnly = TRUE)

v.margins = c(0.1,0.1,0.1,0.1)

nonb = c("all", "APR", "IR", "CRU", "MR", "TRI", "STR", "DR", "SLS",  "GQ", "Z")

l.df.all = {}

####99% Winsorization####

for (k in 1:length(nonb)){
	df.windows = read.csv("results/04_Genome_wide_comparison_to_nonB_motifs/nonB_coverage_versus_read.csv") %>% filter(nonB == nonb[k])
	df.windows$Experiment = paste(df.windows$Species, df.windows$Bio_sample, df.windows$Technique, sep = "_")
	
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
	
	l.df.all[[k]] = bind_rows(l.df2)
}

df.windows = bind_rows(l.df.all)

df.windows$nonB = factor(df.windows$nonB, levels = nonb)

Group = rep("none", nrow(df.windows))
Group[which(grepl("chrX", df.windows$chr))] = "chrX"
Group[which(grepl("chrY", df.windows$chr))] = "chrY"
df.windows$Group = Group

Experiments = unique(df.windows$Experiment)

df = df.windows %>% filter(Experiment %in% c("Ssyndactylus_JAMBI_PDALseq", "Ssyndactylus_Karenina_PDALseq"))
df.IR = df %>% filter(nonB == "IR")
df.GQ = df %>% filter(nonB == "GQ")

df.plot_kmeans = data.frame("chr" = df.IR$chr,"start" = df.IR$start, "end" = df.IR$end, "IR" = df.IR$nonb.count, "GQ" = df.GQ$nonb.count, "PDALseq" = df.IR$reads, "Experiment" = df.IR$Experiment, "interval" = df.IR$interval)

df.Kmeans = df.plot_kmeans  %>% select(IR,GQ,PDALseq) 
KMEANS = kmeans(df.Kmeans, centers = 3, nstart = 20)
df.plot_kmeans$cluster = factor(as.character(KMEANS$cluster), levels = c("1", "2", "3"), labels = c("genome cluster", "alpha-sat cluster", "genome cluster"))


PA = ggplot(df.plot_kmeans, aes(x = IR, y = PDALseq, color = cluster, shape = cluster)) +
	geom_point(alpha = 0.1) +
	theme_classic() +
	facet_wrap(~Experiment, nrow = 1) +
	coord_fixed() +
	xlab("IR density") +
	ylab("PDAL-seq read density") +
	scale_color_manual(values = c("black", "red")) +
	theme(axis.text.x = element_text(angle = 90),
	      legend.position = "bottom")
		

PB = ggplot(df.plot_kmeans, aes(x = GQ, y = PDALseq, color = cluster, shape = cluster)) +
	geom_point(alpha = 0.1) +
	theme_classic() +
	facet_wrap(~Experiment, nrow = 1) +
	coord_fixed() +
	xlab("G4 density") +
	ylab("PDAL-seq read density") +
	scale_color_manual(values = c("black", "red")) +
	theme(axis.text.x = element_text(angle = 90),
	      legend.position = "bottom")
	
PAB = plot_grid(PA, PB, align = "v", nrow = 1, labels = c("A","B"))


#Label in centromere coordinates
df.cen = read.delim("resources/centromeres/Ssyndactylus.bed", header = F)
colnames(df.cen) = c("chr", "start", "end", "annotation")
	
l.df.nonb = df.windows %>% filter(Experiment %in% c("Ssyndactylus_JAMBI_PDALseq")) %>% filter(nonB == "all") %>% group_by(chr) %>% group_split()
          
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
df.sat = read.delim("results/04_Genome_wide_comparison_to_nonB_motifs/Ssyndactylus_subTerm_monhorSF4.bed", header = F)
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
    	ggtitle("subterminal alpha-sat")

P.cluster = ggplot(df.plot_kmeans, aes(x = interval, y = Experiment, fill = cluster)) +
	geom_tile() +
	scale_x_continuous(limits = c(min(df.Sat$interval), max(df.Sat$interval)), expand = c(0, 0)) + scale_fill_manual(values = c("black", "red")) +
	theme_classic() +
    	ggtitle("Cluster from pannels A & B") +
    	theme(legend.position = "bottom")

PC = plot_grid(P.ideogram, P.sat, P.cluster, align = "v", ncol = 1, rel_heights = c(0.8,0.8,3))

P = plot_grid(PAB, PC, ncol = 1, labels = c("", "C"), rel_heights = c(1.2,1))

ggsave(args[1], P, width = 7, height = 4, scale = 1.5)

