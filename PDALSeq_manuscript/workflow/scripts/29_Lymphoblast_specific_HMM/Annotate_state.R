library(tidyverse)
options(scipen = 999)

args       <- commandArgs(trailingOnly = TRUE)
bed_path   <- args[1]
states_csv <- args[2]
output_dir <- args[3]

df.bed   <- readr::read_delim(bed_path, delim = "\t",
                               col_names = c("chr", "start", "end"))[, 1:3]
df.state <- readr::read_csv(states_csv)

df.bed$hidden_state <- df.state$hidden_state

l.df.bed <- df.bed %>% group_by(hidden_state) %>% group_split()

for (i in seq_along(l.df.bed)) {
    STATE <- l.df.bed[[i]]$hidden_state[1]
    PATH  <- file.path(output_dir, paste0("state_", STATE, ".bed"))
    write.table(l.df.bed[[i]], PATH, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
}
