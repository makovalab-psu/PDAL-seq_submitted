library(tidyverse)

state <- c("0", "1", "2", "3", "4", "5", "6", "7")
ann   <- c("Sphase_GM24143", "G0G1_GM24143", "Sphase_GM24385", "G0G1_GM24385")

df <- expand.grid(State = state, ann = ann, stringsAsFactors = FALSE)

Mean.observed <- numeric(nrow(df))
Mean.null     <- numeric(nrow(df))
Wilcox.p      <- numeric(nrow(df))

for (i in seq_len(nrow(df))) {
    path <- paste0(
        "data/25_Sphase_contribution/Multivariate_guassian_HMM/hominoidae/8/",
        "State_functional_element_enrichment/state_", df$State[i], "/",
        df$ann[i], "_enrichment.bed"
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
    "data/25_Sphase_contribution/Multivariate_guassian_HMM/hominoidae/8/State_functional_element_enrichment/Sphase_enrichment.csv",
    row.names = FALSE)
