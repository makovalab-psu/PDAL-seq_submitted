library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

df = readr::read_csv(args[1], show_col_types = FALSE)

PDALseq = c(
	"Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq",
	"Ptroglodytes_AG18354_PDALseq", "Ptroglodytes_AG18358_PDALseq",
	"Ppaniscus_KB8711_PDALseq", "Ppaniscus_KB14571_PDALseq",
	"Ggorilla_KB3781_PDALseq", "Ggorilla_F7549_PDALseq",
	"Ppygmaeus_AG05252_PDALseq", "Pabelii_AG06213_PDALseq",
	"Ssyndactylus_JAMBI_PDALseq", "Ssyndactylus_Karenina_PDALseq",
	"Hsapien_K562_PDALseq"
)

df = df[, which(colnames(df) %in% PDALseq)]
df.clean = df[!rowSums(is.na(df)), ]

write.csv(df.clean, args[2], row.names = FALSE)
