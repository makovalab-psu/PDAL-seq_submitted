library(tidyverse)
library(car)
library(MuMIn)

args <- commandArgs(trailingOnly = TRUE)

options(na.action = "na.fail")

df = readr::read_csv("results/05_nonB_motif_cluster_analysis/nonb_cluster_coverage.csv")

#Normalize the data####

l.df = df %>% group_by(Experiment) %>% group_split()
l.df.norm = {}
for (i in 1:length(l.df)){
	print(i)
	df.i = data.frame(l.df[[i]])
	Norm_reads = (10^9)*df.i$read.count/sum(df.i$read.count)
	df.i$Norm_reads = Norm_reads
	l.df.norm[[i]] = df.i
}

df.norm = bind_rows(l.df.norm)

####inear regression on a down sample###

#Filter the data by cell line
df.i = df.norm %>% filter(Experiment == args[1])
df.i = df.i[sample(1:nrow(df.i), 10000),]

#Fit the data to the global model
fit = lm(Norm_reads ~ length.N + APR + IR + CRU + MR + TRI + Z + STR + DR + SLS + GQ + GC.content, df.i)

#Test all combinations of the parameters	
res <- dredge(fit, trace=2)

#Determine the most probable models
df.sel = data.frame(subset(res, weight >= 0.01, recalc.weights=FALSE))
df.sel$Experiment = args[1]

####Write the model selection table####

write.csv(df.sel, args[2], row.names = TRUE)

####Write the averaged model coefficients####

#define file name
sink(args[3])
summary(model.avg(res))
#close the external connection
sink()

