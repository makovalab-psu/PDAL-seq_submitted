library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

df = read.csv("data/09_Human_HMM/04_Multivariate_guassian/BIC.csv")

P.BIC = ggplot(df, aes(x =States, y = BIC)) +
	geom_point() +
	geom_line() +
	theme_classic() +
	xlab("Hidden States (n)")+
	ylab("Bayesian Information Criterion (BIC)")

ggsave(args[1], P.BIC)

