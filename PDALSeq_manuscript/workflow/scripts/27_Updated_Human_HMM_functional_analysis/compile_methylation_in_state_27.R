library(tidyverse)

state <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")

l.df <- list()

for (i in seq_along(state)) {
    file_path <- paste0(
        "data/27_Updated_Human_HMM_functional_analysis/methylation/state_", state[i], ".txt"
    )
    df.i           <- read.delim(file_path, header = FALSE, sep = "")
    colnames(df.i) <- c("count", "methylation_rate")
    df.i$percent   <- df.i$count / sum(df.i$count)
    df.i$state     <- state[i]
    l.df[[i]]      <- df.i
}

df <- bind_rows(l.df)

write.csv(df,
    "data/27_Updated_Human_HMM_functional_analysis/methylation/hsapiens.csv",
    row.names = FALSE)
