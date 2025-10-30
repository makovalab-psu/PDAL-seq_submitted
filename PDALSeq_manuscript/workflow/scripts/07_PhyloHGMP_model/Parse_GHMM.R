library(tidyverse)
library(viridis)
library(cowplot)

df = readr::read_csv("results/07_PhyloHGMP_model/alignment_clustering/hominoidea/10/GHMM_states.csv")

####Make a tidy formatted data frame####

l.df = {}

for (i in 1:(ncol(df)-1)){
	df.i = data.frame("state" = df[[15]], "reads" = df[[i]], "Experiment" = colnames(df)[i])
	l.df[[i]] = df.i
}

df.tidy = bind_rows(l.df)

####Make a dataframe that summarizes the mean read count in each state####

l.df = {}

for (i in 1:(ncol(df)-1)){
	df.i = data.frame("state" = df[[15]], "reads" = df[[i]])
	l.df[[i]] = data.frame(df.i %>% group_by(state) %>% summarize(MEAN = mean(reads)))
	l.df[[i]]$Experiment = colnames(df)[i] 
}

df.out = bind_rows(l.df)

####Make a data frame to count and arrange the states####

df.arrange = df.tidy %>% group_by(state) %>% summarize(MEAN = mean(reads), COUNT = n()) %>% arrange(MEAN)

####Order the states by mean PDALseq reads####

df.arrange$state = factor(df.arrange$state, levels = df.arrange$state)
df.out$state = factor(df.out$state, levels = df.arrange$state)

####Make a plot that shows the number of blocks in each state####

P.count = ggplot(df.arrange, aes(y = state, x = COUNT)) +
	geom_bar(stat = "identity") +
	theme_classic()

P.Mean = ggplot(df.arrange, aes(x = "", y = state, fill = MEAN)) +
	geom_tile() +
	scale_fill_viridis(option = "turbo") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90))

df.out$Experiment = factor(df.out$Experiment, levels = unique(df.out$Experiment))

P.tile = ggplot(df.out, aes(x = Experiment, y = state, fill = MEAN)) +
	geom_tile() +
	scale_fill_viridis(option = "turbo") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90))

P = plot_grid(P.count, P.Mean, P.tile, nrow = 1, align = "h", rel_widths = c(2,1,5))

