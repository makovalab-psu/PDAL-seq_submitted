library(tidyverse)

args       <- commandArgs(trailingOnly = TRUE)
base_dir   <- args[1]
output_csv <- args[2]

state <- as.character(0:9)
df    <- data.frame(State = state, stringsAsFactors = FALSE)

Mean.observed <- numeric(length(state))
Mean.null     <- numeric(length(state))
Wilcox.p      <- numeric(length(state))

for (i in seq_along(state)) {
    path  <- file.path(base_dir, "State_RNA_expression", paste0("state_", state[i], ".txt"))
    df.i  <- read.csv(path)
    Mean.observed[i] <- mean(df.i$Observed)
    Mean.null[i]     <- mean(df.i$Null)
    Wilcox.p[i]      <- wilcox.test(df.i$Observed, df.i$Null)$p.value
}

df$Mean.observed  <- Mean.observed
df$Mean.null      <- Mean.null
df$Wilcox.p       <- Wilcox.p
df$Enrichment     <- Mean.observed / Mean.null
df$log2Enrichment <- log2(df$Enrichment)

write.csv(df, output_csv, row.names = FALSE)
