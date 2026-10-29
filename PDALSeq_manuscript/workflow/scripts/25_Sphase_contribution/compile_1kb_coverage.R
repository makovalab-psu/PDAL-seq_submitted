library(tidyverse)

dmso_samples <- c(
    "Hsa_GM24143_DMSO_PDALSeq_js4012d20mM_S4_R1_001",
    "Hsa_GM24143_DMSO_PDALSeq_js4012e20mM_S5_R1_001",
    "Hsa_GM24143_DMSO_PDALSeq_js4012f20mM_S6_R1_001",
    "Hsa_GM24385_DMSO_PDALSeq_js4012a20mM_S1_R1_001",
    "Hsa_GM24385_DMSO_PDALSeq_js4012b20mM_S2_R1_001",
    "Hsa_GM24385_DMSO_PDALSeq_js4012c20mM_S3_R1_001"
)

bulk_samples <- c(
    "Hsapien_GM24143_LCL_PDALseq_js2007d20mM_S2_R1_001",
    "Hsapien_GM24143_LCL_PDALseq_js2007e20mM_S3_R1_001",
    "Hsapien_GM24143_LCL_PDALseq_js2007f40mM_S4_R1_001",
    "Hsapien_GM24143_LCL_PDALseq_js2007g40mM_S5_R1_001",
    "Hsapien_GM24385_LCL_PDALseq_js1014A20_1_S2_R1_001",
    "Hsapien_GM24385_LCL_PDALseq_js1014A20_2_S3_R1_001",
    "Hsapien_GM24385_LCL_PDALseq_js1014A40_1_S4_R1_001",
    "Hsapien_GM24385_LCL_PDALseq_js1014A40_2_S5_R1_001"
)

read_bg <- function(path) {
    read.delim(path, sep = "", header = FALSE,
               col.names = c("chr", "start", "end", "reads"))
}

normalize_rpm <- function(df) {
    df$reads / sum(df$reads, na.rm = TRUE) * 1e6
}

coords <- read_bg(paste0("data/25_Sphase_contribution/04_1kb_coverage/dmso/",
                         dmso_samples[1], ".bed"))[, 1:3]

df_out <- coords

for (s in dmso_samples) {
    path <- paste0("data/25_Sphase_contribution/04_1kb_coverage/dmso/", s, ".bed")
    df_out[[s]] <- normalize_rpm(read_bg(path))
}

for (s in bulk_samples) {
    path <- paste0("data/25_Sphase_contribution/04_1kb_coverage/bulk/", s, ".bed")
    df_out[[s]] <- normalize_rpm(read_bg(path))
}

dir.create("results/25_Sphase_contribution", recursive = TRUE, showWarnings = FALSE)
write.csv(df_out, "results/25_Sphase_contribution/Sphase_1kb_coverage.csv",
          row.names = FALSE)
