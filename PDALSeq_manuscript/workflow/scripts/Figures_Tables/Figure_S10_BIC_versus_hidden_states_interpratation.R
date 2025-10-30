library(tidyverse)
library(cowplot)

args = commandArgs(trailingOnly = TRUE)

####BIC####

taxon = c("hominoidae")

l.df = {}
for (i in 1:length(taxon)){
	file_path = paste("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/", taxon[i], "/BIC.csv", sep = "")
	df.i = read.csv(file_path)
	df.i$taxon = taxon[i]
	l.df[[i]] = df.i
}

df = bind_rows(l.df) 

df$taxon = factor(df$taxon, levels = taxon)

P.BIC = ggplot(df, aes(x =States, y = BIC)) +
	geom_point() +
	geom_line() +
	theme_classic() +
	facet_wrap(~taxon, scales = "free") +
	xlab("Hidden States (n)")+
	ylab("Bayesian Information Criterion (BIC)")

####Means####

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
Exp.label   = c("H. sapiens GM24385 PDAL-seq",
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
hidden_states = c("2", "8", "18", "30")


plot_means = function(x = "18"){
	file_path = paste("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/", x, "/Average_PDALseq_in_GHMM_states.csv", sep = "")
	df = readr::read_csv(file_path) %>% arrange(MEAN)
	State.levels = unique(df$state)
	df$state = factor(df$state, levels = State.levels)
	df$Experiment = factor(df$Experiment, levels = Exp.levels, labels = Exp.label)

	P.tile = ggplot(df, aes(x = Experiment, y = state, fill = MEAN)) +
		geom_tile() +
		scale_fill_gradient(low = "white", high = "red") +
		theme_classic() +
		theme(axis.text.x = element_blank(),
	      	      legend.position = "top",
	      	      legend.title = element_blank(),
	      	      axis.text.y = element_blank()) +
		ggtitle(paste("hidden state\nnumber =", x)) +
		ylab("Hidden state")
	
	return(P.tile)
}

l.mean_plots = lapply(hidden_states, plot_means)

####Compile and save plots####

P = plot_grid(l.mean_plots[[1]], l.mean_plots[[2]], l.mean_plots[[3]], l.mean_plots[[4]], ncol = 2)

P.out = plot_grid(P , P.BIC, ncol =1, rel_heights = c(2,1))
ggsave(args[1], P.out, height = 10, width = 6, scale = 1, limitsize = FALSE)
