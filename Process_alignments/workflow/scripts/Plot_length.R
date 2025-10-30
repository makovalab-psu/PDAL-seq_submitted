args = commandArgs(trailingOnly = TRUE)
#args = c("data/check_block_lengths/hominoidae/")

library(tidyverse)

Species = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus")

read.hist = function(x){
    file.path = paste(args[1], x, ".hist", sep = "")
    df = read.delim(file.path, sep = "", header = F)
    colnames(df) = c("Count", "Length")
    df$Species = x
    return(df)
}

df = bind_rows(lapply(Species, read.hist))

df$Species = factor(df$Species, levels = Species, labels = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus"))

P = ggplot(df, aes(x = Length, y = Count, color = Species)) +
    geom_line() +
    theme_classic() +
    scale_x_continuous(trans = "log10")

ggsave(args[2], P, width = 4, height = 3)