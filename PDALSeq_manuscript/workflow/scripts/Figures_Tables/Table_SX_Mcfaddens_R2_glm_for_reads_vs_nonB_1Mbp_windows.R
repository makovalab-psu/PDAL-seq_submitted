library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

df = read.csv(args[1])

df = df %>% select(Experiment, Species, Bio_sample, Technique, mcfadden_r2)

l.df = df %>% group_by(Experiment) %>% group_split()

l.df.2 = {}

for (i in 1:length(l.df)){
	l.df.2[[i]] = data.frame(l.df[[i]][1,])
}

df.2 = bind_rows(l.df.2)

write.csv(df.2, args[2], row.names = F)
