library(tidyverse)
library(viridis)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)
#args = c("data/09_Human_HMM/03_Multivariate_guassian/14/Average_PDALseq_in_GHMM_states.csv", "data/09_Human_HMM/03_Multivariate_guassian/14/Means.csv", "data/09_Human_HMM/03_Multivariate_guassian/14/plot_HMM_vs_real_means.png")

#Read in mean data
df.mean = read.csv(args[1])

#Read in model
df.model = read.csv(args[2]) 

colnames(df.model)[2:(ncol(df.model)-2)] = as.character(unique(df.mean$Experiment))

#Tidy format model data

####Tidy format the df.model data frame####

df = df.model

l.df = {}

for (i in 2:(ncol(df)-2)){
	df.i = data.frame("state" = df[[which(colnames(df) == "state")]], "MEAN" = df[[i]], "Experiment" = colnames(df)[i])
	l.df[[i]] = df.i
}

df.model = bind_rows(l.df)

#### arrange the data ####

df.arrange = df.model %>% group_by(state) %>% summarize(AVERAGE = mean(MEAN)) %>% arrange(AVERAGE)
EXPERIMENT = unique(df.mean$Experiment)
STATES = as.character(df.arrange$state)
df.model$state = factor(df.model$state, levels = STATES)
df.mean$state = factor(df.mean$state, levels = STATES)
df.model$Experiment = factor(df.model$Experiment, levels = EXPERIMENT)
df.mean$Experiment = factor(df.mean$Experiment, levels = EXPERIMENT)

df.model = df.model %>% arrange(state) %>% arrange(Experiment)
df.mean = df.mean %>% arrange(state) %>% arrange(Experiment)


####Make plots####

P.model = ggplot(df.model, aes(x = Experiment, y = state, fill = MEAN)) +
	geom_tile() +
	scale_fill_gradient(low = "white", high = "red") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90),
	      axis.title.x = element_blank(),
	      legend.position = "top") +
	ggtitle("Guassian Mean in state")

P.mean = ggplot(df.mean, aes(x = Experiment, y = state, fill = MEAN)) +
	geom_tile() +
	scale_fill_gradient(low = "white", high = "red") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90),
	      axis.title.x = element_blank(),
	      legend.position = "top") +
	ggtitle("Mean PDAL-seq reads in state")

P.cor = ggplot(mapping = aes(x = df.model$MEAN, y = df.mean$MEAN)) +
	geom_point() +
	coord_fixed() +
	theme_classic() +
	geom_abline(intercept = 0, slope = 1) +
	xlab("Guassian Mean in state") +
	ylab("Mean PDAL-seq reads in state")
	
P = plot_grid(P.model, P.mean, P.cor, nrow = 1)

ggsave(args[3], P, width = 12, height = 5, scale = 1.5)
