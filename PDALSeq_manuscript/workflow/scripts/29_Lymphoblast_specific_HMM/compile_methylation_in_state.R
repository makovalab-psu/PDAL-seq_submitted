library(tidyverse)

args       <- commandArgs(trailingOnly = TRUE)
base_dir   <- args[1]
output_csv <- args[2]

state <- as.character(0:9)
l.df  <- list()

for (i in seq_along(state)) {
    path  <- file.path(base_dir, "methylation", paste0("state_", state[i], ".txt"))
    df.i  <- read.delim(path, header = FALSE, sep = "")
    colnames(df.i) <- c("count", "methylation_rate")
    df.i$percent   <- df.i$count / sum(df.i$count)
    df.i$state     <- state[i]
    l.df[[i]]      <- df.i
}

df <- bind_rows(l.df)
write.csv(df, output_csv, row.names = FALSE)
