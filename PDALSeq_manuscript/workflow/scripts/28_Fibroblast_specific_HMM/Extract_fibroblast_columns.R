library(tidyverse)

args = commandArgs(trailingOnly = TRUE)

df = readr::read_csv(args[1])

PDALseq = c(
    "Ppaniscus_KB8711_PDALseq",
    "Ppaniscus_KB14571_PDALseq",
    "Ggorilla_KB3781_PDALseq",
    "Ggorilla_F7549_PDALseq",
    "Ppygmaeus_AG05252_PDALseq",
    "Pabelii_AG06213_PDALseq"
)

df = df[, which(colnames(df) %in% PDALseq)]
df.clean <- df[!rowSums(is.na(df)), ]
nrow(df.clean)

write.csv(df.clean, args[2], row.names = FALSE)
