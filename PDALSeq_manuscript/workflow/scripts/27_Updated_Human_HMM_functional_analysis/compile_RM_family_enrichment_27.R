library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

State <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
Ann   <- c("Psuedogene", "LINE", "SINE", "LTR", "DNA_transposons", "Satelitte",
           "Simple_repeat", "Low_complexity")

df <- expand.grid(State, Ann)
colnames(df) <- c("State", "Ann")

Mean.observed <- c()
Mean.null     <- c()
Wilcox.p      <- c()

for (i in 1:nrow(df)) {
    File_path <- paste0(
        "data/27_Updated_Human_HMM_functional_analysis/RM_family_enrichment/state_",
        df$State[i], "/", df$Ann[i], "_enrichment.bed"
    )
    df.i <- read.csv(File_path)
    if (length(unique(df.i$Observed)) == 1) {
        Mean.observed[i] <- mean(df.i$Observed)
        Mean.null[i]     <- mean(df.i$Null)
        Test             <- abs(Mean.observed[i] - Mean.null[i])
        Wilcox.p[i]      <- length(which(abs(df.i$Null - Mean.null[i]) >= Test)) / length(df.i$Null)
    } else {
        Mean.observed[i] <- mean(df.i$Observed)
        Mean.null[i]     <- mean(df.i$Null)
        Wilcox.p[i]      <- wilcox.test(df.i$Observed, df.i$Null)$p.value[[1]][1]
    }
}

df$Mean.observed  <- Mean.observed
df$Mean.null      <- Mean.null
df$Wilcox.p       <- Wilcox.p
df$Enrichment     <- Mean.observed / Mean.null
df$log2Enrichment <- log2(df$Enrichment)

write.csv(df, args[1], row.names = FALSE)
