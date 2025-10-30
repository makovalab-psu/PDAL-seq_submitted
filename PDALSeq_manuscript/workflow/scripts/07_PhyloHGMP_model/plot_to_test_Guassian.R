library(tidyverse)
library(cowplot)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)

#args = c("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominidae/20/GHMM_states.csv", "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominidae/20/Means.csv", "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominidae/20/Varience.csv", "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominidae/20/plot_to_test_Guassian.pdf")

df.PDALseq = readr::read_csv(args[1])

N.state = length(unique(df.PDALseq$hidden_state))

df.Mean = read.csv(args[2])
df.Var = read.csv(args[3]) 
Experiments = colnames(df.PDALseq)[1:(ncol(df.PDALseq)-1)]
States = unique(df.PDALseq$hidden_state)

####Determine min/max####

PDALseq = c()

for (i in 1:(ncol(df.PDALseq)-1)){
	PDALseq = c(PDALseq, df.PDALseq[[i]])
}

MIN = min(PDALseq)
MAX = max(PDALseq)

####Make histograms####

l.row = {}


for (i in 1:length(States)){
	State = States[i]
	df.i = df.PDALseq %>% filter(hidden_state == State)
	df.Var.i = df.Var %>% filter(state == State)
	df.Mean.i = df.Mean %>% filter(state == State)

	Experiment = Experiments[i]

	l.Hist = {}
	Mu = c()
	Var = c()

	for (j in 1:(ncol(df.PDALseq)-1)){
		Mu[j] = df.Mean.i[[j+1]]
		print(Mu[j])
		Var[j] = df.Var.i[[j+1]]
		print(Var[j])
		n <- nrow(df.i)
		bin_width <- diff(range(df.i[[j]])) / 30  # Assuming 30 bins
		l.P.hist = {}
		Hist = ggplot() +
			geom_histogram(mapping = aes(x = df.i[[j]], y = ..density..), fill = "dimgrey", color = "white") +
			stat_function(fun = dnorm, args = list(mean = Mu[j], sd = sqrt(Var[j])), color = "red", size = 1) +
			theme_classic() +
			labs(title = paste(colnames(df.i)[j], "\nstate =", State), x = "PDAL-seq", y = "Count")
		l.P.hist[[1]] = Hist
		l.Hist[[j]] = grid.arrange(grobs = l.P.hist, ncol = 1)
	}
	

	l.row[[i]] = grid.arrange(grobs = l.Hist, ncol = length(l.Hist))

}


P.out = grid.arrange(grobs = l.row, ncol = 1)

ggsave(args[4], plot = P.out, width = 14, height = 1.5*N.state, scale = 5, limitsize = FALSE)

