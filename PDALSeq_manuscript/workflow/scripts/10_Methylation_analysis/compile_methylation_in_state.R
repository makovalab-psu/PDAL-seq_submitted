library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

Species = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii")
State = c("0", "1", "2", "3", "4", "5", "6", "7")

df = data.frame(expand.grid(Species, State))
colnames(df) = c("Species", "State")

Start = c(0,0.2,0.4,0.6,0.8)
End = c(0.2,0.4,0.6,0.8,1.0)

#Start = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
#End = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)


l.df = {}
for (i in 1:nrow(df)){
	file_path = paste("data/10_Methylation_analysis/02_state_annotations/hominoidae/8/", df$Species[i], "/state_", df$State[i], ".txt", sep = "")
	paste(file_path)
	df.i = read.delim(file_path, head = F, col.names = c("Count", "CpG_rate"), sep = "")
	Count = c()
	for (j in 1:length(Start)){
		if (End[j] == 1.0){
			Count[j] = sum(df.i$Count[which(df.i$CpG_rate >= Start[j] & df.i$CpG_rate <= End[j])])
		}else{
			Count[j] = sum(df.i$Count[which(df.i$CpG_rate >= Start[j] & df.i$CpG_rate < End[j])])
		}
	}
	Sum = sum(Count)
	Percent = 100*Count/Sum
	Species = df$Species[i]
	State = df$State[i]
	l.df[[i]] = data.frame(Species, State, "CpG_5mC_rate_start" = Start, "CpG_5mC_rate_end" =  End, Count, Percent)
	print(file_path)
	print(Sum == sum(df.i$Count))
}

df.out = bind_rows(l.df)
df.out$Range = paste(df.out$CpG_5mC_rate_start, df.out$CpG_5mC_rate_end, sep = "_")

write.csv(df.out, args[1], row.names = F) 
