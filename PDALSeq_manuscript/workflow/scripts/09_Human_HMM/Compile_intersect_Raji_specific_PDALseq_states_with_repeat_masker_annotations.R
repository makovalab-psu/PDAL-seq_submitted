library(dplyr)

Class = c("DNA","DNA?","DNA/Crypton","DNA/Crypton-A","DNA/hAT","DNA/hAT?","DNA/hAT-Ac","DNA/hAT-Blackjack","DNA/hAT-Charlie","DNA/hAT-hAT19","DNA/hAT-Tag1","DNA/hAT-Tip100","DNA/hAT-Tip100?","DNA?/hAT-Tip100?","DNA/Kolobok","DNA/Merlin","DNA/MULE-MuDR","DNA/PIF-Harbinger","DNA/PiggyBac","DNA/TcMar","DNA/TcMar-Mariner","DNA/TcMar-Pogo","DNA/TcMar-Tc1","DNA/TcMar-Tc2","DNA/TcMar-Tigger","LINE/CR1","LINE/Dong-R4","LINE/I-Jockey","LINE/L1","LINE/L1-Tx1","LINE/L2","LINE/Penelope","LINE/RTE-BovB","LINE/RTE-X","Low_complexity","LTR","LTR/ERV1","LTR/ERVK","LTR/ERVL","LTR/ERVL-MaLR","LTR/Gypsy","RC/Helitron","Retroposon/SVA","rRNA","Satellite","Satellite/acro","Satellite/centr","Satellite/subtelo","scRNA","Simple_repeat","SINE/5S-Deu-L2","SINE/Alu","SINE/MIR","SINE/tRNA","SINE/tRNA-Deu","SINE/tRNA-RTE","snRNA","srpRNA","tRNA","Unknown")

State_files = c("State_7_all", "State_7_high_PDAL-seq_expected", "State_7_variable_PDAL-seq_expected", "State_7_low_PDAL-seq_expected", "State_7_nonsyntenic", "State_7_syntenic")

df = data.frame(expand.grid(Class, State_files))
colnames(df) = c("Class", "State")

Mean.observed = c()
Mean.null = c()
Wilcox.p = c()
Parent_class = c()
Subclass = c()

l.Class = strsplit(as.character(df$Class), split = "/") 


for (i in 1:nrow(df)){
	File_path = paste("data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/repeat_masker_intersects/", df$State[i], "/", gsub("/", "-", df$Class[i]), ".csv", sep = "")
	df.i = read.csv(File_path)
	
	if (length(unique(df.i$Observed)) <= 1){
		Mean.observed[i] = df.i$Observed[1]
		Mean.null[i] = mean(df.i$Null)
		Xnull = abs(df.i$Null - Mean.null[i])
		Xobs = abs(Mean.observed[i] - Mean.null[i])
		Wilcox.p[i] = length(which(Xnull >= Xobs))/length(Xnull)
	}else{
		Mean.observed[i] = mean(df.i$Observed)
		Mean.null[i] = mean(df.i$Null)
		Wilcox = wilcox.test(df.i$Observed, df.i$Null)
		Wilcox.p[i] = Wilcox$p.value[[1]][1]
	}
	print(paste(df$State[i], df$Class[i], Mean.observed[i], Mean.null[i],  Wilcox.p[i]))
	Parent_class[i] = l.Class[[i]][1]
	if (length(l.Class[[i]]) > 1){Subclass[i] = l.Class[[i]][2]}else{Subclass[i] = "none"}
}


df$Mean.observed = Mean.observed
df$Mean.null = Mean.null
df$Enrichment = Mean.observed/Mean.null
df$log2Enrichment = log2(Mean.observed/Mean.null)
df$Wilcox.p = Wilcox.p
df$Parent_class = Parent_class
df$Subclass = Subclass

write.csv(df, "data/09_Human_HMM/05_analysis_of_specific_states/Human_state_7/repeat_masker_intersects.csv", row.names = FALSE)
