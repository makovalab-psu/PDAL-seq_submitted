library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
out_file <- args[1]

# Match compile_unique_map_estimation.R logic

genome_size <- function(chrom_sizes_path) {
  df <- read.delim(chrom_sizes_path, header = FALSE)
  sum(df$V2)
}

bed_size <- function(bed_path) {
  df <- read.delim(bed_path, header = FALSE)
  sum(df$V3 - df$V2 + 1)
}

calc_umr <- function(path) {
  df     <- read.delim(path, header = TRUE, sep = "")
  hq_idx <- which(df$Primary.score.minus.suplemental.score >= 10)
  uniquely.mapped <- if (length(hq_idx) > 0) sum(df$count[hq_idx]) else sum(df$count)
  uniquely.mapped / sum(df$count)
}

# All new experiments use the Hsapien genome/blacklist
total_len     <- genome_size("resources/genomes/Hsapien.chrom.sizes")
blacklist_len <- bed_size("results/03_Correlation_and_black_list_unmappable_regions/blacklist_files/Hsapien.bed")

BLACKLIST_EXP <- c(
  "Hsapien_K562_PDALseq",
  "Hsapien_GM24143_DMSO_PDALseq",
  "Hsapien_GM24385_DMSO_PDALseq"
)

exp_meta <- tribble(
  ~Experiment,                       ~Species,  ~Bio_sample, ~Technique,
  "Hsapien_K562_PDALseq",            "Hsapien", "K562",      "PDALseq",
  "Hsapien_GM24143_DMSO_PDALseq",    "Hsapien", "GM24143",   "PDALseq",
  "Hsapien_GM24385_DMSO_PDALseq",    "Hsapien", "GM24385",   "PDALseq",
)

base <- "data/30_update_SI_tables_with_revisions_data/Blacklist_unique_maps"

rows <- lapply(BLACKLIST_EXP, function(exp) {
  inside  <- calc_umr(file.path(base, "inside",  paste0(exp, ".txt")))
  outside <- calc_umr(file.path(base, "outside", paste0(exp, ".txt")))
  meta    <- exp_meta[exp_meta$Experiment == exp, ]
  data.frame(
    Species                    = meta$Species,
    Bio_sample                 = meta$Bio_sample,
    Technique                  = meta$Technique,
    Experiment                 = exp,
    Blacklist.length           = blacklist_len,
    Fraction_genome_in_blacklist = blacklist_len / total_len,
    Unique.maps.blacklist        = inside,
    Unique.maps.outside.blacklist = outside
  )
})

df <- bind_rows(rows)

SPECIES       <- c("Hsapien")
SPECIES.label <- c("H. sapiens")
df$Species <- factor(df$Species, levels = SPECIES, labels = SPECIES.label)

df <- df %>%
  arrange(Bio_sample) %>%
  arrange(Species) %>%
  arrange(Technique)

write.csv(df, out_file, row.names = FALSE)
