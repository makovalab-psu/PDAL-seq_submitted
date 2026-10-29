library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
out_file <- args[1]

# Helpers — match Compile_mapping_statistics_by_biolgical_sample.R logic

calc_flag <- function(path) {
  df <- read.delim(path, header = FALSE, sep = "")
  TOTAL      <- sum(df$V1[df$V2 %in% c(0L, 4L, 16L)])
  unmapped   <- ifelse(4L %in% df$V2, df$V1[df$V2 == 4L], 0L)
  mapped     <- sum(df$V1[df$V2 %in% c(0L, 16L)])
  secondary  <- sum(df$V1[df$V2 %in% c(2048L, 2064L)])
  data.frame(TOTAL = TOTAL, unmapped = unmapped, mapped = mapped, secondary = secondary)
}

calc_qmap <- function(path) {
  df <- read.delim(path, header = FALSE, sep = "")
  data.frame(
    less.than.equal.to.Q10 = sum(df$V1[df$V2 <= 10]),
    greater.than.Q10       = sum(df$V1[df$V2 >= 10]),
    greater.than.Q20       = sum(df$V1[df$V2 >= 20]),
    greater.than.Q30       = sum(df$V1[df$V2 >= 30])
  )
}

calc_umr <- function(path) {
  df     <- read.delim(path, header = TRUE, sep = "")
  hq_idx <- which(df$Primary.score.minus.suplemental.score >= 10)
  uniquely.mapped <- if (length(hq_idx) > 0) sum(df$count[hq_idx]) else sum(df$count)
  data.frame(uniquely.mapped = uniquely.mapped)
}

# Experiment metadata (hardcoded — K562/DMSO are not in PDALseq_samples.txt)
K562_EXP <- c("Hsapien_K562_PDALseq_0mM", "Hsapien_K562_PDALseq_20mM",
               "Hsapien_K562_PDALseq_40mM", "Hsapien_K562_PDALseq")
DMSO_EXP <- c("Hsapien_GM24143_DMSO_PDALseq", "Hsapien_GM24385_DMSO_PDALseq")

exp_meta <- tribble(
  ~Experiment,                         ~Species,   ~Bio_sample, ~Technique,
  "Hsapien_K562_PDALseq_0mM",          "Hsapien",  "K562",      "PDALseq",
  "Hsapien_K562_PDALseq_20mM",         "Hsapien",  "K562",      "PDALseq",
  "Hsapien_K562_PDALseq_40mM",         "Hsapien",  "K562",      "PDALseq",
  "Hsapien_K562_PDALseq",              "Hsapien",  "K562",      "PDALseq",
  "Hsapien_GM24143_DMSO_PDALseq",      "Hsapien",  "GM24143",   "PDALseq",
  "Hsapien_GM24385_DMSO_PDALseq",      "Hsapien",  "GM24385",   "PDALseq",
)

# K562: All_reads + High_MapQ; DMSO: All_reads only (no MapQ filter in step 25)
exp_mapq <- bind_rows(
  data.frame(Experiment = K562_EXP, MapQ = "All_reads",  stringsAsFactors = FALSE),
  data.frame(Experiment = K562_EXP, MapQ = "High_MapQ",  stringsAsFactors = FALSE),
  data.frame(Experiment = DMSO_EXP, MapQ = "All_reads",  stringsAsFactors = FALSE)
)

base <- "data/30_update_SI_tables_with_revisions_data"

rows <- lapply(seq_len(nrow(exp_mapq)), function(i) {
  exp <- exp_mapq$Experiment[i]
  mq  <- exp_mapq$MapQ[i]
  flag <- calc_flag(file.path(base, "Merged_bam_flags",           mq, paste0(exp, ".txt")))
  qual <- calc_qmap(file.path(base, "Merged_mapping_quality",     mq, paste0(exp, ".txt")))
  umr  <- calc_umr( file.path(base, "Merged_primary_vs_supplemental", mq, paste0(exp, ".txt")))
  meta <- exp_meta[exp_meta$Experiment == exp, ]
  data.frame(
    Species          = meta$Species,
    Bio_sample       = meta$Bio_sample,
    Technique        = meta$Technique,
    Experiment       = exp,
    MapQ             = mq,
    mapped           = flag$mapped,
    secondary        = flag$secondary,
    Total            = flag$TOTAL,
    uniquely.mapped  = umr$uniquely.mapped,
    Fraction.unique  = umr$uniquely.mapped / flag$TOTAL,
    less.than.equal.to.Q10 = qual$less.than.equal.to.Q10,
    greater.than.Q10       = qual$greater.than.Q10,
    greater.than.Q20       = qual$greater.than.Q20,
    greater.than.Q30       = qual$greater.than.Q30
  )
})

df <- bind_rows(rows)

SPECIES       <- c("Hsapien")
SPECIES.label <- c("H. sapiens")
df$Species <- factor(df$Species, levels = SPECIES, labels = SPECIES.label)
df$MapQ    <- factor(df$MapQ, levels = c("All_reads", "High_MapQ"),
                     labels = c("All reads", "High MapQ"))

df <- df %>%
  arrange(Bio_sample) %>%
  arrange(Species) %>%
  arrange(Technique) %>%
  arrange(MapQ)

write.csv(df, out_file, row.names = FALSE)
