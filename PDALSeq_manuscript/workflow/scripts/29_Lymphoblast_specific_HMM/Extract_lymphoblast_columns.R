library(tidyverse)

args = commandArgs(trailingOnly = TRUE)

df = readr::read_csv(args[1])

PDALseq = c(
    "Hsapien_GM24385_PDALseq",
    "Hsapien_GM24143_PDALseq",
    "Ptroglodytes_AG18354_PDALseq",
    "Ptroglodytes_AG18358_PDALseq",
    "Ssyndactylus_JAMBI_PDALseq",
    "Ssyndactylus_Karenina_PDALseq"
)

df = df[, which(colnames(df) %in% PDALseq)]
df.clean <- df[!rowSums(is.na(df)), ]
nrow(df.clean)

write.csv(df.clean, args[2], row.names = FALSE)
