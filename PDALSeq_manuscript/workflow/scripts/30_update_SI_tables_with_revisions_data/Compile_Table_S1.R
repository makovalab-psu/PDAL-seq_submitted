library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
out_file <- args[1]

# Hardcoded sample metadata — new samples are not in PDALseq_samples.txt.
# Trimmed_fastq and Mapped_bam (pre-dedup) columns are NA: those intermediate
# files are temp() in steps 21/25 and are deleted after their downstream rules
# complete. Only Raw_fastq, Mapped_deduplicated_bam, and (K562 only)
# Mapped_deduplicated_high_MapQ_bam are available without re-running alignment.
sample_meta <- tribble(
  ~Sequencing_sample,                                   ~Species,   ~Bio_sample, ~Technique, ~Experiment,
  "Hsa_K562_LCL_PDALSeq_js4003_0mM_S1_R1_001",         "Hsapien",  "K562",      "PDALseq",  "Hsapien_K562_PDALseq_0mM",
  "Hsa_K562_LCL_PDALSeq_js4003_20mM1_S2_R1_001",       "Hsapien",  "K562",      "PDALseq",  "Hsapien_K562_PDALseq_20mM",
  "Hsa_K562_LCL_PDALSeq_js4003_20mM2_S3_R1_001",       "Hsapien",  "K562",      "PDALseq",  "Hsapien_K562_PDALseq_20mM",
  "Hsa_K562_LCL_PDALSeq_js4003_40mM1_S4_R1_001",       "Hsapien",  "K562",      "PDALseq",  "Hsapien_K562_PDALseq_40mM",
  "Hsa_K562_LCL_PDALSeq_js4003_40mM2_S5_R1_001",       "Hsapien",  "K562",      "PDALseq",  "Hsapien_K562_PDALseq_40mM",
  "Hsa_GM24143_DMSO_PDALSeq_js4012d20mM_S4_R1_001",    "Hsapien",  "GM24143",   "PDALseq",  "Hsapien_GM24143_DMSO_PDALseq",
  "Hsa_GM24143_DMSO_PDALSeq_js4012e20mM_S5_R1_001",    "Hsapien",  "GM24143",   "PDALseq",  "Hsapien_GM24143_DMSO_PDALseq",
  "Hsa_GM24143_DMSO_PDALSeq_js4012f20mM_S6_R1_001",    "Hsapien",  "GM24143",   "PDALseq",  "Hsapien_GM24143_DMSO_PDALseq",
  "Hsa_GM24385_DMSO_PDALSeq_js4012a20mM_S1_R1_001",    "Hsapien",  "GM24385",   "PDALseq",  "Hsapien_GM24385_DMSO_PDALseq",
  "Hsa_GM24385_DMSO_PDALSeq_js4012b20mM_S2_R1_001",    "Hsapien",  "GM24385",   "PDALseq",  "Hsapien_GM24385_DMSO_PDALseq",
  "Hsa_GM24385_DMSO_PDALSeq_js4012c20mM_S3_R1_001",    "Hsapien",  "GM24385",   "PDALseq",  "Hsapien_GM24385_DMSO_PDALseq",
)

K562_SAMPLES <- c(
  "Hsa_K562_LCL_PDALSeq_js4003_0mM_S1_R1_001",
  "Hsa_K562_LCL_PDALSeq_js4003_20mM1_S2_R1_001",
  "Hsa_K562_LCL_PDALSeq_js4003_20mM2_S3_R1_001",
  "Hsa_K562_LCL_PDALSeq_js4003_40mM1_S4_R1_001",
  "Hsa_K562_LCL_PDALSeq_js4003_40mM2_S5_R1_001"
)

read_count <- function(path) {
  as.integer(readLines(path, n = 1))
}

read_flags <- function(path) {
  if (!file.exists(path)) return(data.frame(mapped = NA_integer_, supplementary = NA_integer_))
  df <- read.delim(path, header = FALSE, sep = "")
  data.frame(
    mapped        = sum(df$V1[df$V2 %in% c(0L, 16L)]),
    supplementary = sum(df$V1[df$V2 %in% c(2048L, 2064L)])
  )
}

rows <- lapply(seq_len(nrow(sample_meta)), function(i) {
  s    <- sample_meta$Sequencing_sample[i]
  raw  <- read_count(paste0("data/30_update_SI_tables_with_revisions_data/Raw_fastq/", s, ".txt"))
  dedup <- read_flags(paste0("data/30_update_SI_tables_with_revisions_data/Mapped_deduplicated_bam/", s, ".txt"))
  if (s %in% K562_SAMPLES) {
    hq <- read_flags(paste0("data/30_update_SI_tables_with_revisions_data/Mapped_deduplicated_high_MapQ_bam/", s, ".txt"))
  } else {
    hq <- data.frame(mapped = NA_integer_, supplementary = NA_integer_)
  }
  data.frame(
    Species                                             = sample_meta$Species[i],
    Bio_sample                                          = sample_meta$Bio_sample[i],
    Technique                                           = sample_meta$Technique[i],
    Sequencing_sample                                   = s,
    Experiment                                          = sample_meta$Experiment[i],
    Raw.fastq.million.reads                             = raw,
    Trimmed.fastq.million.reads                         = NA_integer_,
    Mapped.bam.primary.million.reads                    = NA_integer_,
    Mapped.bam.supplementary.million.reads              = NA_integer_,
    Mapped.deduplicated.bam.primary.million.reads       = dedup$mapped,
    Mapped.deduplicated.bam.supplementary.million.reads = dedup$supplementary,
    Mapped.deduplicated.high.MapQ.bam.primary.million.reads       = hq$mapped,
    Mapped.deduplicated.high.MapQ.bam.supplementary.million.reads = hq$supplementary
  )
})

df <- bind_rows(rows)

SPECIES       <- c("Hsapien")
SPECIES.label <- c("H. sapiens")
df$Species <- factor(df$Species, levels = SPECIES, labels = SPECIES.label)

df <- df %>%
  arrange(Bio_sample) %>%
  arrange(Technique)

write.csv(df, out_file, row.names = FALSE)
