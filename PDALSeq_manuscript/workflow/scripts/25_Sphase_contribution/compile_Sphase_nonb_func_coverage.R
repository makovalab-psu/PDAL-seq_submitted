library(tidyverse)

ann_list <- c("Sphase_GM24143", "G0G1_GM24143", "Sphase_GM24385", "G0G1_GM24385")
feat_list <- c(
    "promoter", "utr5", "exon", "intron", "utr3",
    "censat", "repeats",
    "APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z",
    "coreorigins", "stochorigins"
)

enrich_base <- "data/25_Sphase_contribution/06_nonb_func_enrichment"
ref_base    <- "data/25_Sphase_contribution/05_merged_windows"

df <- expand.grid(ann = ann_list, feat = feat_list, stringsAsFactors = FALSE)

df$bp_overlap <- mapply(function(an, fe) {
    path <- file.path(enrich_base, an, paste0(fe, "_bp_overlap.txt"))
    as.numeric(readLines(path))
}, df$ann, df$feat)

ref_total_bp <- setNames(
    sapply(ann_list, function(an) {
        bed <- read.table(file.path(ref_base, paste0(an, ".bed")), sep = "\t")
        sum(bed[, 3] - bed[, 2])
    }),
    ann_list
)

df$ref_total_bp     <- ref_total_bp[df$ann]
df$fraction_covered <- df$bp_overlap / df$ref_total_bp

write.csv(df,
    "results/25_Sphase_contribution/Sphase_nonb_func_coverage.csv",
    row.names = FALSE)
