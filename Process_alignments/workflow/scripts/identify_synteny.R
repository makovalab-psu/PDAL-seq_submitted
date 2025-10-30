
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)

#args = c("data/10_quantify_syteny/Hsapien/chr4.txt", "0.99", "data/11_identify_syteny/Hsapien/chr4.txt")
Threshold = as.numeric(args[2])

df = read.delim(args[1])
df$species = factor(df$species, levels = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus"), labels = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus"))

l.df = {}
Species = unique(df$species)

for (i in 1:length(Species)){
    df.i = df %>% filter(species == Species[i]) %>% group_by(chr)
    df.i = data.frame(table(df.i$chr)) %>% arrange(-Freq)
    df.i$Var1 = factor(df.i$Var1, levels = df.i$Var1)
    df.i$species = Species[i]
    df.i$N = 1:nrow(df.i)
    cum.percent = c()
    Have.not.hit.threshold = T
    Keep = c()
    for (k in 1:nrow(df.i)){
        cum.percent[k] = sum(df.i$Freq[1:k])/sum(df.i$Freq)
        Keep[k] = Have.not.hit.threshold
        if (cum.percent[k] >= Threshold){Have.not.hit.threshold = F}
    }
    df.i$cum.percent = cum.percent
    df.i$Keep = Keep
    l.df[[i]] = df.i
}

df.cumulative = bind_rows(l.df)


P.synteny = ggplot(df, aes(x = block, y = start, color = chr)) +
    geom_point(alpha = 0.1) +
    theme_classic() +
    facet_wrap(~species) +
    theme(legend.position = "none")

P.cumulative = ggplot(df.cumulative, aes(x = N, y = cum.percent, color = Keep)) +
    geom_point(alpha = 1) +
    theme_classic() +
    facet_wrap(~species) +
    theme() +
    xlab("Chromosome index") +
    ylab("Cumulative percent") +
    scale_y_continuous(limits = c(0, 1))

P.count = ggplot(df.cumulative, aes(x = N, y = Freq, color = Keep)) +
    geom_point(alpha = 1) +
    theme_classic() +
    facet_wrap(~species) +
    theme() +
    xlab("Chromosome index") +
    ylab("Repesentatives in alignment percent")

df.out = df.cumulative %>% filter(Keep) %>% select(species, Var1)

write.table(df.out, args[3], col.names = F, row.names = F, sep = "\t", quote = F)

ggsave(paste(args[3], "synteny.png", sep = "_"), P.synteny, width = 12, height = 10)
ggsave(paste(args[3], "sumulative.png", sep = "_"), P.cumulative, width = 12, height = 10)
ggsave(paste(args[3], "count.png", sep = "_"), P.count, width = 12, height = 10)