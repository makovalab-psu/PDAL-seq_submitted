library(tidyverse)

state <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")

df <- data.frame(State = state)

Mean.observed <- c()
Mean.null     <- c()
Wilcox.p      <- c()

for (i in 1:nrow(df)) {
    File_path <- paste0(
        "data/27_Updated_Human_HMM_functional_analysis/State_RNA_expression/state_",
        df$State[i], ".txt"
    )
    df.i          <- read.csv(File_path)
    Mean.observed[i] <- mean(df.i$Observed)
    Mean.null[i]     <- mean(df.i$Null)
    Wilcox.p[i]      <- wilcox.test(df.i$Observed, df.i$Null)$p.value[[1]][1]
}

df$Mean.observed  <- Mean.observed
df$Mean.null      <- Mean.null
df$Wilcox.p       <- Wilcox.p
df$Enrichment     <- Mean.observed / Mean.null
df$log2Enrichment <- log2(df$Enrichment)

write.csv(df,
    "data/27_Updated_Human_HMM_functional_analysis/State_RNA_expression/hsapiens.csv",
    row.names = FALSE)
