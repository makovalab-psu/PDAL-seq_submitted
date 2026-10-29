library(tidyverse)

existing_experiments = c(
	"Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq",
	"Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq",
	"Ptroglodytes_AG18354_PDALseq", "Ptroglodytes_AG18358_PDALseq",
	"Ppaniscus_KB8711_PDALseq", "Ppaniscus_KB14571_PDALseq",
	"Ggorilla_KB3781_PDALseq", "Ggorilla_F7549_PDALseq",
	"Ppygmaeus_AG05252_PDALseq", "Pabelii_AG06213_PDALseq",
	"Ssyndactylus_JAMBI_PDALseq", "Ssyndactylus_Karenina_PDALseq",
	"Hsapien_GM24385_TruSeq", "Ptroglodytes_AG18354_TruSeq",
	"Ppaniscus_KB8711_TruSeq", "Ggorilla_KB3781_TruSeq",
	"Ppygmaeus_AG05252_TruSeq", "Pabelii_AG06213_TruSeq",
	"Ssyndactylus_JAMBI_TruSeq"
)

l.M = list()

for (i in seq_along(existing_experiments)){
	csv_i = paste0("data/07_PhyloHGMP_model/02_Normalize_alignment_coverage/hominoidae/",
	               existing_experiments[i], ".csv")
	df.i = readr::read_csv(csv_i, show_col_types = FALSE)
	M.i = matrix(df.i$Normalized_read_count_per_bp_per_billion_reads, nrow = 1)
	rownames(M.i) = existing_experiments[i]
	l.M[[i]] = M.i
}

df.k562 = readr::read_csv(
	"data/26_Ape_PCA_and_Human_HMM_with_K562/02_normalize_alignment_coverage/hominoidae/Hsapien_K562_PDALseq.csv",
	show_col_types = FALSE)
M.k562 = matrix(df.k562$Normalized_read_count_per_bp_per_billion_reads, nrow = 1)
rownames(M.k562) = "Hsapien_K562_PDALseq"
l.M[[length(l.M) + 1]] = M.k562

M.out = do.call(rbind, l.M)

dir.create("data/26_Ape_PCA_and_Human_HMM_with_K562/compile_alignment_coverage/hominoidae",
           recursive = TRUE, showWarnings = FALSE)
write.csv(t(M.out),
          "data/26_Ape_PCA_and_Human_HMM_with_K562/compile_alignment_coverage/hominoidae/alignment_coverage_matrix.csv",
          row.names = FALSE)
