library(dplyr)
library(ggplot2)
library(viridis)

args <- commandArgs(trailingOnly = TRUE)

df = readr::read_delim(args[1])


States = unique(df$hidden_state)
States = States[order(States)]

Experiments = colnames(df)[1:(ncol(df)-1)]
Sample_1 = c()
Sample_2 = c()
for (i in 1:(length(Experiments)-1)){
	y = Experiments[(i+1):length(Experiments)]
	x = rep(Experiments[i], length(y))
	Sample_1 = c(Sample_1, x)
	Sample_2 = c(Sample_2, y)
}

df.index = data.frame(Sample_1,Sample_2)

l.df = {}
#Mean_spearman = c()
Mean_euclidean_distance = c()
Mean_conservation = c()

for (i in 1:length(States)){
	print(States[i])
	df.i = df.index
	df.state = df %>% filter(hidden_state == States[i])
	#Spearman = c()
	Euclidean_distance = c()
	nrow(df.i)
	for (j in 1:nrow(df.i)){
		print(j)
		print(paste(df.i$Sample_1[j], df.i$Sample_2[j]))
		x = df.state[[which(colnames(df.state) == df.i$Sample_1[j])]]
		y = df.state[[which(colnames(df.state) == df.i$Sample_2[j])]]
		x = x/mean(x)
		y = y/mean(y)
		#test.j = cor.test(x, y, method = "spearman")
		#Spearman[j] = test.j$estimate
		Euclidean_distance[j] <- sqrt(sum((x - y)^2))
		print(Euclidean_distance[j])
		#print(Spearman[j])
		
	}
	#df.i$Spearman = Spearman
	df.i$Euclidean_distance = Euclidean_distance
	df.i$Conservation = 1/(1+Euclidean_distance)
	#df.i$Mean_spearman[i] = prod(Spearman)^(1/length(Spearman))
	Mean_euclidean_distance[i] = mean(df.i$Euclidean_distance)
	Mean_conservation[i] = mean(df.i$Conservation)
	df.i$State = States[i]
	l.df[[i]] = df.i
}

df.out = data.frame(States, Mean_euclidean_distance, Mean_conservation)

df.individual_state = bind_rows(l.df)
df.individual_state$Sample_1 = factor(df.individual_state$Sample_1, levels = Experiments)
df.individual_state$Sample_2 = factor(df.individual_state$Sample_2, levels = Experiments)

P = ggplot(df.individual_state, aes(x = Sample_1, y = Sample_2, fill = Conservation)) + 
	geom_tile() +
	facet_wrap(~State) +
	scale_fill_viridis(trans = "log10") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))

write.csv(df.out, args[2], row.names = F)
write.csv(df.individual_state , args[3], row.names = F)
ggsave(args[4], P, width = 11, height = 11)
	
