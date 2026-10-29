#!/usr/bin/env Rscript
# compile_results.R
# Compile per-region read counts into a single CSV with metadata.
#
# Usage:
#   Rscript compile_results.R <output_csv>
#
# Output columns:
#   Dataset         experiment name (e.g. Hsapien_GM24385_PDALseq)
#   Total_bases     total aligned reference bases across primary mapped reads
#   Syntenic        "Syntenic" or "Non-syntenic"
#   Region          feature category (promoter, utr5, exon, intron, utr3,
#                   censat, nongenic_repeats, nongenic_nonrepetitive)
#   Region_nt       nucleotides in the region BED
#   Bases_in_region read-aligned bases falling inside the region
#   Species         species name
#   Method          PDALseq or TruSeq (parsed from experiment name)

args <- commandArgs(trailingOnly = TRUE)
output_file <- args[1]

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# ---------------------------------------------------------------------------
# Configuration (must match 19_PDAL_Seq_signal_in_and_out_of_syntenic_regions.smk)
# ---------------------------------------------------------------------------
syntenic_species <- c("Ggorilla", "Hsapien", "Ppaniscus", "Ptroglodytes", "Pabelii", "Ppygmaeus")

regions <- c(
  "syntenic_promoter", "syntenic_utr5", "syntenic_exon",
  "syntenic_intron",   "syntenic_utr3", "syntenic_censat",
  "syntenic_nongenic_repeats", "syntenic_nongenic_nonrepetitive",
  "nonsyntenic_promoter", "nonsyntenic_utr5", "nonsyntenic_exon",
  "nonsyntenic_intron",   "nonsyntenic_utr3", "nonsyntenic_censat",
  "nonsyntenic_nongenic_repeats", "nonsyntenic_nongenic_nonrepetitive"
)

# ---------------------------------------------------------------------------
# Load sample metadata
# ---------------------------------------------------------------------------
samples <- read_tsv("PDALseq_samples.txt", col_types = cols(.default = "c")) %>%
  select(Experiment, Species) %>%
  distinct()

# Derive Method from experiment name (last underscore-delimited token)
samples <- samples %>%
  mutate(Method = str_extract(Experiment, "[^_]+$"))

# Keep only the 6 target species
samples <- samples %>% filter(Species %in% syntenic_species)

# ---------------------------------------------------------------------------
# Read all count files
# ---------------------------------------------------------------------------
results <- list()

for (sp in syntenic_species) {
  sp_exps <- samples %>% filter(Species == sp) %>% pull(Experiment)

  for (exp in sp_exps) {
    for (region in regions) {
      path <- file.path(
        "data/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/02_read_counts",
        sp, region, paste0(exp, ".txt")
      )
      if (!file.exists(path)) {
        warning("Missing: ", path)
        next
      }
      vals <- read_csv(path, col_names = c("Total_bases", "Bases_in_region", "Region_nt"),
                       col_types = "ddd", show_col_types = FALSE)
      vals$Dataset  <- exp
      vals$Species  <- sp
      vals$Region_full <- region
      results[[length(results) + 1]] <- vals
    }
  }
}

result_df <- bind_rows(results)

# ---------------------------------------------------------------------------
# Parse Syntenic / Region columns from Region_full
# ---------------------------------------------------------------------------
result_df <- result_df %>%
  mutate(
    Syntenic = if_else(str_starts(Region_full, "syntenic_"), "Syntenic", "Non-syntenic"),
    Region   = str_remove(Region_full, "^(syntenic|nonsyntenic)_")
  ) %>%
  select(-Region_full)

# Join method from samples metadata
result_df <- result_df %>%
  left_join(samples %>% select(Experiment, Method), by = c("Dataset" = "Experiment"))

# ---------------------------------------------------------------------------
# Final column order matching the requested output
# ---------------------------------------------------------------------------
result_df <- result_df %>%
  select(Dataset, Total_bases, Syntenic, Region, Region_nt, Bases_in_region, Species, Method)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_csv(result_df, output_file)
message("Written: ", output_file)
