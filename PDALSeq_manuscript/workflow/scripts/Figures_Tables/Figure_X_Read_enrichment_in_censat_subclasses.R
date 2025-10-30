library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)
library(ggbeeswarm)

CenSat_levels = c("active_hor", "hor", "dhor", "mon-hor", "cenSat", "ct", "bSat", "gSat", "HSat2", "HSat2_3", "HSat3", "subTerm", "mixedSuperFamily")

nonb_levels = c("all", "APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS", "GQ")

Exp.levels   = c("Hsapien_GM24385_PDALseq",
		 "Hsapien_GM24385_TruSeq",
		 "Ptroglodytes_AG18354_PDALseq",
		 "Ptroglodytes_AG18354_TruSeq",
		 "Ppaniscus_KB8711_PDALseq",
		 "Ppaniscus_KB8711_TruSeq",
		 "Ggorilla_KB3781_PDALseq",
		 "Ggorilla_KB3781_TruSeq",
		 "Ppygmaeus_AG05252_PDALseq",
		 "Ppygmaeus_AG05252_TruSeq",
		 "Pabelii_AG06213_PDALseq",
		 "Pabelii_AG06213_TruSeq",
		 "Ssyndactylus_JAMBI_PDALseq",
		 "Ssyndactylus_JAMBI_TruSeq")

Exp.labels   = c("H. sapiens GM24385 PDAL-seq",
		 "H. sapiens GM24385 TruSeq",
		 "P. troglodytes AG18354 PDAL-seq",
		 "P. troglodytes AG18354 TruSeq",
		 "P. paniscus KB8711 PDAL-seq",
		 "P. paniscus KB8711 TruSeq",
		 "G. gorilla KB3781 PDAL-seq",
		 "G. gorilla KB3781 TruSeq",
		 "P. pygmaeus AG05252 PDAL-seq",
		 "P. pygmaeus AG05252 TruSeq",
		 "P. abelii AG06213 PDAL-seq",
		 "P. abelii AG06213 TruSeq",
		 "S. syndactylus JAMBI PDAL-seq",
		 "S. syndactylus JAMBI TruSeq")
		 
Species_levels = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Pabelii", "Ppygmaeus", "Ssyndactylus")	
Species_labels = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. abelii", "P. pygmaeus", "S. syndactylus")		 

df.enr = read.delim("results/17_PDALseq_reads_in_CenSat/Censat_subclass_enrichment.txt", sep = "\t") %>% filter(!is.na(Enrichment)) %>% arrange(-Enrichment)


hist(nchar(unique(df.enr$Subclass)))
Long_names = unique(df.enr$Subclass)[which(nchar(unique(df.enr$Subclass)) > 30)]
New_names = paste("Subclass", 1:length(Long_names))
for (i in 1:length(Long_names)){
	df.enr$Subclass[which(df.enr$Subclass == Long_names[i])] = New_names[i]
}
hist(nchar(unique(df.enr$Subclass)))
hist(df.enr$Enrichment)

df.enr$Enrichment[which(df.enr$Enrichment > 3)] = 3
MIN_max = c(0, 3)

l.df.enr = df.enr %>% group_by(Class) %>% group_split()
l.p = {}
rows = c()

hist(nchar(unique(df.i$Subclass)))
for (i in 1:length(l.df.enr)){
	df.i = l.df.enr[[i]] %>% arrange(-Enrichment)
	rows[i] = length(unique(df.i$Subclass))
	df.i$ylab = factor(df.i$Subclass, levels = rev(unique(df.i$Subclass)))
	df.i$Experiment = factor(df.i$Experiment, levels = Exp.levels, labels = Exp.labels)
	Sig = rep("ns", nrow(df.i))
	Sig[which(df.i$Wilcox.p <= 0.05)] = "*"
	Sig[which(df.i$Wilcox.p <= 0.01)] = "**"

	Enr_labels = paste(round(df.i$Enrichment, digits = 1), Sig, sep = "\n")
	Enr_labels[which(is.na(df.i$Enrichment))] = "ns"
	df.i$Enr_labels = Enr_labels
	l.p[[i]] = ggplot(df.i, aes(x = Experiment, y = ylab, fill = Enrichment, label = Enr_labels)) +
		ggtitle(df.i$Class[1]) +
		geom_tile() +
		geom_text(color = "white", size = 3) +
		scale_fill_viridis(option = "turbo", limits = MIN_max) +
		theme_classic() +
		theme(axis.text.x =  element_text(angle = 45, hjust = 1),
	      		axis.text.y =  element_text(angle = 45, vjust = 1),
	      		legend.position = "right",
	      		axis.title.x = element_blank(),
	      		axis.title.y = element_blank()) +
	      	scale_x_discrete(limits = Exp.labels)
	
}

rows + 3

P.grid = plot_grid(l.p[[1]], l.p[[2]], l.p[[3]], l.p[[4]], l.p[[5]], l.p[[6]], l.p[[7]], l.p[[8]], l.p[[9]], l.p[[10]], l.p[[11]], l.p[[12]], ncol = 1, align = "v", rel_heights = rows+4)


ggsave("plots/Figure_X_Read_enrichment_in_censat_subclasses/Figure_X.svg", P.grid, width = 6.5, height = 0.5*sum(rows +3), units = "in", limitsize = FALSE)

write.table(data.frame(Long_names, New_names), "plots/Figure_X_Read_enrichment_in_censat_subclasses/Long_names.txt", quote = F)

