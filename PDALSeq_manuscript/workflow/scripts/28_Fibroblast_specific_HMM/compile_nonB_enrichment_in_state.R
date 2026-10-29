library(tidyverse)

args       <- commandArgs(trailingOnly = TRUE)
base_dir   <- args[1]
output_csv <- args[2]

state <- as.character(0:9)
nonb  <- c("APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z")

df <- expand.grid(State = state, nonB = nonb, stringsAsFactors = FALSE)

Mean.observed <- numeric(nrow(df))
Mean.null     <- numeric(nrow(df))
Wilcox.p      <- numeric(nrow(df))

for (i in seq_len(nrow(df))) {
    path <- file.path(base_dir, "State_nonB_enrichment",
                      paste0("state_", df$State[i]),
                      paste0(df$nonB[i], "_enrichment.bed"))
    df.i          <- read.csv(path)
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
