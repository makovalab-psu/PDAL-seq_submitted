library(tidyverse)

args       <- commandArgs(trailingOnly = TRUE)
base_dir   <- args[1]
output_csv <- args[2]

state <- as.character(0:9)
ann   <- c("cds_protCode_singleTrans", "ngnr_protCode_singleTrans", "chrG_cpgislands",
           "promoter_protCode_singleTrans", "enhancer_protCode_singleTrans",
           "repeats_protCode_singleTrans", "exon_protCode_singleTrans",
           "rna_protCode_singleTrans", "gene_nonprotCode", "utr3_protCode_singleTrans",
           "gene_protCode_singleTrans", "utr5_protCode_singleTrans",
           "intron_protCode_singleTrans", "Recombination",
           "core_stochastic_origins_chm13v2.0_hglft.ucsc",
           "core_origins_hglft_genome_e0596_1cc9e0",
           "stochastic_origins_hglft_genome_2ff9ed_1ccfd0")

df <- expand.grid(State = state, ann = ann, stringsAsFactors = FALSE)

Mean.observed <- numeric(nrow(df))
Mean.null     <- numeric(nrow(df))
Wilcox.p      <- numeric(nrow(df))

for (i in seq_len(nrow(df))) {
    path <- file.path(base_dir, "State_functional_element_enrichment",
                      paste0("state_", df$State[i]),
                      paste0(df$ann[i], "_enrichment.bed"))
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
