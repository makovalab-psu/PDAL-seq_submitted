library(tidyverse)

####A plot of PDAL-seq reads versus hsa locus####

df.chr = readr::read_csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/GHMM_states.csv")
df.bg = readr::read_delim("data/07_PhyloHGMP_model/01_alignment_coverage/hominoidae/Hsapien_GM24385_PDALseq.bg", col_names = c("chr", "start", "end", "reads"))[,1:3]
df.G4 = readr::read_delim("data/04_Genome_wide_comparison_to_nonB_motifs/03_nonb_density/10000_nucleotides/Hsapien/GQ.bg", delim = "\t", col_names = c("chr", "start", "end", "count", "BP", "Coverage", "density"))
df.TRI = readr::read_delim("data/04_Genome_wide_comparison_to_nonB_motifs/03_nonb_density/10000_nucleotides/Hsapien/TRI.bg", delim = "\t", col_names = c("chr", "start", "end", "count", "BP", "Coverage", "density"))

df.chr = cbind(df.chr, df.bg)

range(df.chr[,1:12])


df.chr$state = factor(df.chr$hidden_state, levels = State.levels)

plot_hmm = function(CHR = "chr19", xlim = c(-1.3,5), ylim = c(1.5*10^6, 2*10^6)){
	df.i = df.chr %>% filter(chr == CHR)
	df.i$state = factor(as.character(df.i$hidden_state), levels = State.levels)
	df.G4.i = df.G4 %>% filter(chr == CHR)
	df.TRI.i = df.TRI %>% filter(chr == CHR)
	XLIMITS = xlim
	YLIMITS = ylim
	
	P.G4 = ggplot(df.G4.i, aes(ymin = start, ymax = end, xmin = 0, xmax =  1, fill = density)) +
		geom_rect() +
		theme_classic()  +
		scale_fill_gradient(low = "white", high = "red") +
		scale_y_continuous(limits = YLIMITS, expand = c(0,0)) +
		ggtitle("G4")+
		#theme_void() +
		theme(legend.position = "none")
	P.TRI = ggplot(df.TRI.i, aes(ymin = start, ymax = end, xmin = 0, xmax =  1, fill = density)) +
		geom_rect() +
		theme_classic()  +
		scale_fill_gradient(low = "white", high = "red") +
		scale_y_continuous(limits = YLIMITS, expand = c(0,0)) +
		ggtitle("TRI")+
		#theme_void() +
		theme(legend.position = "none")
		
	P.human = ggplot(df.i, aes(ymin = start, ymax = end, xmin = -1.1, xmax =  Hsapien_GM24385_PDALseq, fill = state)) +
		geom_rect() +
		theme_classic()  +
		scale_x_continuous(limits = XLIMITS, expand = c(0,0)) +
		scale_y_continuous(limits = YLIMITS, expand = c(0,0)) +
		ggtitle("human") +
		theme(legend.position = "left")
		
	plot_grid(P.human, P.G4, P.TRI, align = "v", nrow = 1, rel_widths = c(1.2, 0.5, 0.5))
	
	P.chimp = ggplot(df.i, aes(xmin = start, xmax = end, ymin = -1, ymax =  Ptroglodytes_AG18354_PDALseq, fill = state)) +
		geom_rect() +
		theme_classic()  +
		scale_x_continuous(limits = XLIMITS, expand = c(0,0)) +
		scale_y_continuous(limits = YLIMITS, expand = c(0,0)) +
		ggtitle("chimp")
	P.bonobo = ggplot(df.i, aes(xmin = start, xmax = end, ymin = -1, ymax =  Ppaniscus_KB14571_PDALseq, fill = state)) +
		geom_rect() +
		theme_classic()  +
		scale_x_continuous(limits = XLIMITS, expand = c(0,0)) +
		scale_y_continuous(limits = YLIMITS, expand = c(0,0)) +
		ggtitle("bonobo")
	P.bonobo = ggplot(df.i, aes(xmin = start, xmax = end, ymin = -1, ymax =  Ppaniscus_KB14571_PDALseq, fill = state)) +
		geom_rect() +
		theme_classic()  +
		scale_x_continuous(limits = XLIMITS, expand = c(0,0)) +
		scale_y_continuous(limits = YLIMITS, expand = c(0,0)) +
		ggtitle("bonobo")
	P.gorilla = ggplot(df.i, aes(xmin = start, xmax = end, ymin = -1, ymax =  Ggorilla_F7549_PDALseq, fill = state)) +
		geom_rect() +
		theme_classic()  +
		scale_x_continuous(limits = XLIMITS, expand = c(0,0)) +
		scale_y_continuous(limits = YLIMITS, expand = c(0,0)) +
		ggtitle("gorilla")
	P.orangutan = ggplot(df.i, aes(xmin = start, xmax = end, ymin = -1, ymax =  Ppygmaeus_AG05252_PDALseq, fill = state)) +
		geom_rect() +
		theme_classic()  +
		scale_x_continuous(limits = XLIMITS, expand = c(0,0)) +
		scale_y_continuous(limits = YLIMITS, expand = c(0,0)) +
		ggtitle("orangutan")
	P.siamang = ggplot(df.i, aes(xmin = start, xmax = end, ymin = -1, ymax =  Ssyndactylus_JAMBI_PDALseq, fill = state)) +
		geom_rect() +
		theme_classic()  +
		scale_x_continuous(limits = XLIMITS, expand = c(0,0)) +
		scale_y_continuous(limits = YLIMITS, expand = c(0,0)) +
		ggtitle("siamang")
	P.out = plot_grid(P.human, P.chimp, P.bonobo, P.gorilla, P.orangutan, P.siamang, ncol = 1, align = "v")
	return(P.out)
}


