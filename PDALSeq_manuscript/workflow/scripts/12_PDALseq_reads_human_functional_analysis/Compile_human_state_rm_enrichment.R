
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

State = c("state_0", "state_1", "state_2", "state_3", "state_4", "state_5", "state_6", "state_7", "state_8", "state_9")
Ann = c("DNA","DNA?","DNA/Crypton","DNA/Crypton-A","DNA/hAT","DNA/hAT?","DNA/hAT-Ac","DNA/hAT-Blackjack","DNA/hAT-Charlie","DNA/hAT-hAT19","DNA/hAT-Tag1","DNA/hAT-Tip100","DNA/hAT-Tip100?","DNA?/hAT-Tip100?","DNA/Kolobok","DNA/Merlin","DNA/MULE-MuDR","DNA/PIF-Harbinger","DNA/PiggyBac","DNA/TcMar","DNA/TcMar-Mariner","DNA/TcMar-Pogo","DNA/TcMar-Tc1","DNA/TcMar-Tc2","DNA/TcMar-Tigger","LINE/CR1","LINE/Dong-R4","LINE/I-Jockey","LINE/L1","LINE/L1-Tx1","LINE/L2","LINE/Penelope","LINE/RTE-BovB","LINE/RTE-X","Low_complexity","LTR","LTR/ERV1","LTR/ERVK","LTR/ERVL","LTR/ERVL-MaLR","LTR/Gypsy","RC/Helitron","Retroposon/SVA","rRNA","Satellite","Satellite/acro","Satellite/centr","Satellite/subtelo","scRNA","Simple_repeat","SINE/5S-Deu-L2","SINE/Alu","SINE/MIR","SINE/tRNA","SINE/tRNA-Deu","SINE/tRNA-RTE","snRNA","srpRNA","tRNA","Unknown")

df = expand.grid(State, Ann)
colnames(df) = c("State", "Ann")

Mean.observed = c()
Mean.null = c()
Wilcox.p = c()

for (i in 1:nrow(df)){
	File_path = paste("data/12_PDALseq/human_state_rm_enrichment/", df$Ann[i], "/", df$State[i], ".bed", sep = "")
	df.i = read.csv(File_path)
	if (length(unique(df.i$Observed)) == 1){
		Mean.observed[i] = mean(df.i$Observed)
		Mean.null[i] = mean(df.i$Null)
		Test = abs(Mean.observed[i] - Mean.null[i])
		P = length(which(abs(df.i$Null - Mean.null[i]) >= Test))/length(df.i$Null)
		Wilcox.p[i] = P
	}else{
		Mean.observed[i] = mean(df.i$Observed)
		Mean.null[i] = mean(df.i$Null)
		Wilcox = wilcox.test(df.i$Observed, df.i$Null)
		Wilcox.p[i] = Wilcox$p.value[[1]][1]
	}
}

df$Mean.observed = Mean.observed
df$Mean.null = Mean.null 
df$Wilcox.p = Wilcox.p
df$Enrichment = Mean.observed/Mean.null
df$log2Enrichment = log2(df$Enrichment)

write.csv(df, args[1], row.names = F)
