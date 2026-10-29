#!/usr/bin/env Rscript
# compile_results.R
# Compile per-experiment strand-specific coverage profiles into a single CSV.
#
# Usage:
#   Rscript compile_results.R <output_csv>
#
# Reads all profile TSVs from:
#   data/20_Strand_specific_analysis/03_profiles/{species}/{nonb}/{experiment}.tsv
#   data/20_Strand_specific_analysis/03_profiles/{species}/{nonb}/{experiment}_shuffled.tsv
#
# Output columns:
#   Dataset              experiment name
#   Species              species name
#   Method               PDALseq or TruSeq (last underscore-delimited token)
#   NonB_type            non-B DNA motif type (APR, GQ, etc.)
#   Position             bin center in bp relative to feature 5' end (-995 to 995)
#   Forward_coverage     mean + strand read CPM; NA for strand-aware (BED6) features
#   Reverse_coverage     mean − strand read CPM; NA for strand-aware (BED6) features
#   Template_coverage    mean template-strand CPM; NA for unstranded (BED3) features
#   Nontemplate_coverage mean nontemplate-strand CPM; NA for unstranded (BED3) features
#   Strand_aware         TRUE when template/nontemplate columns are populated (GQ only);
#                        FALSE when forward/reverse columns are populated (all other non-B types)
#
# Strand_aware = TRUE (GQ, BED6):
#   Template_coverage    = mean CPM of displaced ssDNA strand reads near the G4 motif
#   Nontemplate_coverage = mean CPM of G4-forming strand reads
#   Forward_coverage / Reverse_coverage = NA
#
# Strand_aware = FALSE (all other non-B types, BED3):
#   Forward_coverage = mean CPM of + genomic strand reads
#   Reverse_coverage = mean CPM of − genomic strand reads
#   Template_coverage / Nontemplate_coverage = NA

args <- commandArgs(trailingOnly = TRUE)
output_file <- args[1]

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
nonb_types   <- c("APR", "CRU", "DR", "GQ", "TRI", "Z")
profile_base <- "data/20_Strand_specific_analysis/03_profiles"

# ---------------------------------------------------------------------------
# Load sample metadata
# ---------------------------------------------------------------------------
experiments_20 <- c(
  "Hsapien_GM24385_PDALseq",      "Hsapien_GM24385_TruSeq",
  "Ggorilla_KB3781_PDALseq",      "Ggorilla_KB3781_TruSeq",
  "Ptroglodytes_AG18354_PDALseq", "Ptroglodytes_AG18354_TruSeq",
  "Ppaniscus_KB8711_PDALseq",     "Ppaniscus_KB8711_TruSeq",
  "Ppygmaeus_AG05252_PDALseq",    "Ppygmaeus_AG05252_TruSeq",
  "Pabelii_AG06213_PDALseq",      "Pabelii_AG06213_TruSeq",
  "Ssyndactylus_JAMBI_PDALseq",   "Ssyndactylus_JAMBI_TruSeq"
)

samples <- read_tsv("PDALseq_samples.txt", col_types = cols(.default = "c")) %>%
  select(Experiment, Species) %>%
  distinct() %>%
  filter(Experiment %in% experiments_20) %>%
  mutate(Method = str_extract(Experiment, "[^_]+$"))

species_list <- unique(samples$Species)

# ---------------------------------------------------------------------------
# Parse a single profile TSV.
# Handles two formats:
#   BED3 (unstranded):    columns "forward", "reverse"
#   BED6 (strand-aware):  columns "template", "nontemplate"
# Returns a data.frame with Position plus the relevant coverage columns and
# a Strand_aware flag. Returns NULL if the file is empty or unparseable.
# ---------------------------------------------------------------------------
parse_profile <- function(path) {
  if (!file.exists(path) || file.size(path) == 0) return(NULL)

  lines <- readLines(path)
  data_lines <- lines[!startsWith(lines, "#")]
  if (length(data_lines) < 2) return(NULL)

  tbl <- read_tsv(
    paste(data_lines, collapse = "\n"),
    col_types = cols(.default = "c"),
    show_col_types = FALSE
  )

  nms     <- names(tbl)
  pos_col <- nms[1]

  # Detect profile type. Use anchored regex so "nontemplate" does not match
  # the "template" pattern.
  t_col  <- nms[str_detect(nms, regex("^template$",    ignore_case = TRUE))]
  nt_col <- nms[str_detect(nms, regex("^nontemplate$", ignore_case = TRUE))]
  fwd_col <- nms[str_detect(nms, regex("forward",      ignore_case = TRUE))]
  rev_col <- nms[str_detect(nms, regex("reverse",      ignore_case = TRUE))]

  strand_aware <- length(t_col) > 0 && length(nt_col) > 0

  if (strand_aware) {
    tbl %>%
      select(
        Position             = all_of(pos_col),
        Template_coverage    = all_of(t_col[1]),
        Nontemplate_coverage = all_of(nt_col[1])
      ) %>%
      mutate(
        Position             = as.numeric(Position),
        Template_coverage    = as.numeric(Template_coverage),
        Nontemplate_coverage = as.numeric(Nontemplate_coverage),
        Forward_coverage     = NA_real_,
        Reverse_coverage     = NA_real_,
        Strand_aware         = TRUE
      ) %>%
      filter(!is.na(Position))
  } else if (length(fwd_col) > 0 && length(rev_col) > 0) {
    tbl %>%
      select(
        Position         = all_of(pos_col),
        Forward_coverage = all_of(fwd_col[1]),
        Reverse_coverage = all_of(rev_col[1])
      ) %>%
      mutate(
        Position             = as.numeric(Position),
        Forward_coverage     = as.numeric(Forward_coverage),
        Reverse_coverage     = as.numeric(Reverse_coverage),
        Template_coverage    = NA_real_,
        Nontemplate_coverage = NA_real_,
        Strand_aware         = FALSE
      ) %>%
      filter(!is.na(Position))
  } else {
    NULL
  }
}

# ---------------------------------------------------------------------------
# Iterate over all species × nonB × experiments
# ---------------------------------------------------------------------------
results <- list()

for (sp in species_list) {
  sp_exps <- samples %>% filter(Species == sp) %>% pull(Experiment)

  for (nonb in nonb_types) {
    for (exp in sp_exps) {
      path <- file.path(profile_base, sp, nonb, paste0(exp, ".tsv"))

      prof <- parse_profile(path)
      if (is.null(prof)) {
        warning("Missing or empty: ", path)
        next
      }

      prof$Dataset   <- exp
      prof$Species   <- sp
      prof$NonB_type <- nonb
      results[[length(results) + 1]] <- prof
    }
  }
}

result_df <- bind_rows(results)

result_df <- result_df %>%
  left_join(samples %>% select(Experiment, Method), by = c("Dataset" = "Experiment")) %>%
  select(Dataset, Species, Method, NonB_type, Position,
         Forward_coverage, Reverse_coverage,
         Template_coverage, Nontemplate_coverage,
         Strand_aware)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_csv(result_df, output_file)
message("Written: ", output_file)
