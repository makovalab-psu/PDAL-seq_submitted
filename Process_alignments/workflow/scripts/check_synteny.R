
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)

#args = c("data/13_quantify_syteny_after_enforcement/Hsapien/chr2.txt")

df = read.delim(args[1])
df$species = factor(df$species, levels = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus"), labels = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus"))
df$chr.label = paste(df$species, df$chr)

P.synteny = ggplot(df, aes(x = block, y = start, color = chr.label)) +
    geom_point() +
    theme_classic() +
    facet_wrap(~species) 


ggsave(args[2], P.synteny, width = 12, height = 10)