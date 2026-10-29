library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

State <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
Ann   <- c("DNA","DNA?","DNA/Crypton","DNA/Crypton-A","DNA/hAT","DNA/hAT?","DNA/hAT-Ac",
           "DNA/hAT-Blackjack","DNA/hAT-Charlie","DNA/hAT-hAT19","DNA/hAT-Tag1",
           "DNA/hAT-Tip100","DNA/hAT-Tip100?","DNA?/hAT-Tip100?","DNA/Kolobok",
           "DNA/Merlin","DNA/MULE-MuDR","DNA/PIF-Harbinger","DNA/PiggyBac","DNA/TcMar",
           "DNA/TcMar-Mariner","DNA/TcMar-Pogo","DNA/TcMar-Tc1","DNA/TcMar-Tc2",
           "DNA/TcMar-Tigger","LINE/CR1","LINE/Dong-R4","LINE/I-Jockey","LINE/L1",
           "LINE/L1-Tx1","LINE/L2","LINE/Penelope","LINE/RTE-BovB","LINE/RTE-X",
           "Low_complexity","LTR","LTR/ERV1","LTR/ERVK","LTR/ERVL","LTR/ERVL-MaLR",
           "LTR/Gypsy","RC/Helitron","Retroposon/SVA","rRNA","Satellite","Satellite/acro",
           "Satellite/centr","Satellite/subtelo","scRNA","Simple_repeat","SINE/5S-Deu-L2",
           "SINE/Alu","SINE/MIR","SINE/tRNA","SINE/tRNA-Deu","SINE/tRNA-RTE","snRNA",
           "srpRNA","tRNA","Unknown")

df <- expand.grid(State, Ann)
colnames(df) <- c("State", "Ann")

l.df <- list()

for (i in 1:nrow(df)) {
    File_path <- paste0(
        "data/27_Updated_Human_HMM_functional_analysis/human_state_rm_intersect_bp/",
        df$Ann[i], "/", df$State[i], ".bed"
    )
    df.i       <- read.csv(File_path)
    df.i$Ann   <- df$Ann[i]
    df.i$State <- paste0("state_", df$State[i])
    l.df[[i]]  <- df.i
}

df.out <- bind_rows(l.df)
write.csv(df.out, args[1], row.names = FALSE)
