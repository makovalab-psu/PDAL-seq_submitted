library(tidyverse)

ann_list <- c("Sphase_GM24143", "G0G1_GM24143", "Sphase_GM24385", "G0G1_GM24385")
feat_list <- c(
    "promoter", "utr5", "exon", "intron", "utr3",
    "censat", "repeats",
    "APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z",
    "coreorigins", "stochorigins"
)

df <- expand.grid(ann = ann_list, feat = feat_list, stringsAsFactors = FALSE)

Mean.observed <- numeric(nrow(df))
Mean.null     <- numeric(nrow(df))
Wilcox.p      <- numeric(nrow(df))

for (i in seq_len(nrow(df))) {
    path <- paste0(
        "data/25_Sphase_contribution/06_nonb_func_enrichment/",
        df$ann[i], "/", df$feat[i], "_enrichment.bed"
    )
    d <- read.csv(path)
    Mean.observed[i] <- mean(d$Observed)
    Mean.null[i]     <- mean(d$Null)
    Wilcox.p[i]      <- wilcox.test(d$Observed, d$Null)$p.value
}

df$Mean.observed  <- Mean.observed
df$Mean.null      <- Mean.null
df$Wilcox.p       <- Wilcox.p
df$Enrichment     <- Mean.observed / Mean.null
df$log2Enrichment <- log2(df$Enrichment)

write.csv(df,
    "results/25_Sphase_contribution/Sphase_nonb_func_enrichment.csv",
    row.names = FALSE)
