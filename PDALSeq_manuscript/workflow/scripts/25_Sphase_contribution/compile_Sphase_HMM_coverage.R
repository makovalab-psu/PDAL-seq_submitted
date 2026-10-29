library(tidyverse)

state_list <- c("0", "1", "2", "3", "4", "5", "6", "7")
ann_list   <- c("Sphase_GM24143", "G0G1_GM24143", "Sphase_GM24385", "G0G1_GM24385")

enrich_base <- "data/25_Sphase_contribution/Multivariate_guassian_HMM/hominoidae/8/State_functional_element_enrichment"
state_base  <- "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8"

df <- expand.grid(State = state_list, ann = ann_list, stringsAsFactors = FALSE)

df$bp_overlap <- mapply(function(st, an) {
    path <- file.path(enrich_base, paste0("state_", st), paste0(an, "_bp_overlap.txt"))
    as.numeric(readLines(path))
}, df$State, df$ann)

state_total_bp <- setNames(
    sapply(state_list, function(st) {
        bed <- read.table(file.path(state_base, paste0("state_", st, ".bed")), sep = "\t")
        sum(bed[, 3] - bed[, 2])
    }),
    state_list
)

df$state_total_bp   <- state_total_bp[df$State]
df$fraction_covered <- df$bp_overlap / df$state_total_bp

write.csv(df, file.path(enrich_base, "Sphase_HMM_coverage.csv"), row.names = FALSE)
