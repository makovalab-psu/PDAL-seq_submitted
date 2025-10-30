args = commandArgs(trailingOnly = TRUE)

library(tidyverse)

meadianL = function(x = "Hsapien"){
    file_path = paste(args[1], x, ".hist", sep = "")
    df = read.delim(file_path, header = F, sep = "")
    colnames(df) = c("count", "length")
    df = df %>% arrange(length)
    cum_sum = rep(NA, nrow(df))
    for (i in 1:nrow(df)){
        cum_sum[i] = sum(df$count[1:i])
    }
    median = df$length[which.min(abs(cum_sum - 0.5*sum(df$count)))]
    return(median)
}

Species = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus")

MedianL = rep(NA, length(Species))

for (i in 1:length(Species)){
    MedianL[i] = meadianL(Species[i])
}

df = data.frame(Species, MedianL)

write.table(df, args[2], row.names = F, quote = F)