options(scipen=999)

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

df = read_delim(args[1], delim = "\t", col_names = F, show_col_types = FALSE)
colnames(df) = c("chr", "start", "end", "nonb", "cluster")
l.df = df %>% group_by(cluster) %>% group_split()

Lines = rep("", length(l.df))

pb <- txtProgressBar(min = 1,      # Minimum value of the progress bar
                     max = length(l.df), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     char = "=")

for (i in 1:length(l.df)){
    #setTxtProgressBar(pb, i)
    chr = l.df[[i]]$chr[1]
    start = l.df[[i]]$start[1]
    end = l.df[[i]]$end[nrow(l.df[[i]])]
    cluster.ID = l.df[[i]]$cluster[1]
    length.N = end - start +1
    MR = length(which(l.df[[i]]$nonb == "MR"))
    TRI = length(which(l.df[[i]]$nonb == "TRI"))
    STR = length(which(l.df[[i]]$nonb == "STR"))
    Z = length(which(l.df[[i]]$nonb == "Z"))
    APR = length(which(l.df[[i]]$nonb == "APR"))
    DR = length(which(l.df[[i]]$nonb == "DR"))
    SLS = length(which(l.df[[i]]$nonb == "SLS"))
    GQ = length(which(l.df[[i]]$nonb == "GQ"))
    IR = length(which(l.df[[i]]$nonb == "IR"))
    CRU = length(which(l.df[[i]]$nonb == "CRU"))
    Lines[i] = paste(chr, start, end, cluster.ID, length.N, MR, TRI, STR, Z, APR, DR, SLS, GQ, IR, CRU, sep = "\t")
}

CON = file(args[2])
writeLines(Lines, CON)
close(CON)
