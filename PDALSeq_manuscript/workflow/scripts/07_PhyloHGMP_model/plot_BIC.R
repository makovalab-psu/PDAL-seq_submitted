library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

#args = c("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominini/BIC.csv", "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominini/BIC_plot.svg")

df = read.csv(args[1])

P = ggplot(df, aes(x = States, y = BIC)) +
	geom_point() +
	geom_line() +
	theme_classic() +
	xlab("number of hidden states") +
	ylab("BIC")
	
ggsave(args[2], P, width = 5, height = 4, units = "in", scale = 1)
