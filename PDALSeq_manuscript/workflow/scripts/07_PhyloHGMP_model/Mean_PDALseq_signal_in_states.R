library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

df = readr::read_csv(args[1])

####Make a dataframe that summarizes the mean read count in each state####

l.df = {}

for (i in 1:(ncol(df)-1)){
	df.i = data.frame("state" = df[[ncol(df)]], "reads" = df[[i]])
	l.df[[i]] = data.frame(df.i %>% group_by(state) %>% summarize(MEAN = mean(reads)))
	l.df[[i]]$Experiment = colnames(df)[i] 
}

df.out = bind_rows(l.df)

write.csv(df.out, args[2], row.names = F)
