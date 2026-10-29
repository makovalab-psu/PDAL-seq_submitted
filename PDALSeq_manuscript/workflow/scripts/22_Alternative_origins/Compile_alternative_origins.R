library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

all_experiments <- c(
    "Hsapien_CCL86-Raji_PDALseq", "Hsapien_GM24385_PDALseq",
    "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq",
    "Hsapien_K562_PDALseq"
)
shuffle <- c("unscrambled", "scrambled")
state   <- c("all", "4")
origin  <- c(
    "GSM5658908_Ini-seq2.called.replication.origins",
    "GSM5658909_SNS.seq.fractions.0.5.2.kb_called.replication.origins",
    "GSM5658910_SNS.seq.fractions.2.4.kb_called.replication.origins",
    "GSM6127777_Ini-seq2.called.replication.origin.15.min"
)

df <- expand.grid(all_experiments, shuffle, state, origin, stringsAsFactors = FALSE)
colnames(df) <- c("sample", "shuffle", "state", "origin")

l.df <- vector("list", nrow(df))
for (i in seq_len(nrow(df))) {
    file_path <- paste0(
        "data/22_Alternative_origins/Mean_signal/", df$sample[i],
        "/", df$shuffle[i],
        "/state_", df$state[i],
        "/", df$origin[i], ".csv"
    )
    df.i          <- read.csv(file_path)
    df.i$sample   <- df$sample[i]
    df.i$shuffle  <- df$shuffle[i]
    df.i$state    <- df$state[i]
    df.i$origin   <- df$origin[i]
    l.df[[i]]     <- df.i
}

df.out <- bind_rows(l.df)
write.csv(df.out, args[1], row.names = FALSE)
