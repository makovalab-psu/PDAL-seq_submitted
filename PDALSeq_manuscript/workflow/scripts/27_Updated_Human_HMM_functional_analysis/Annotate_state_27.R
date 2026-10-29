library(tidyverse)

options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)

df.bed   <- readr::read_delim(args[1], delim = "\t", col_names = c("chr", "start", "end"))[, 1:3]
df.state <- readr::read_csv(args[2])

df.bed$hidden_state <- df.state$hidden_state

l.df.bed <- df.bed %>% group_by(hidden_state) %>% group_split()

for (i in seq_along(l.df.bed)) {
    STATE <- l.df.bed[[i]]$hidden_state[1]
    PATH  <- paste0("data/27_Updated_Human_HMM_functional_analysis/state_", STATE, ".bed")
    write.table(l.df.bed[[i]], PATH, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
