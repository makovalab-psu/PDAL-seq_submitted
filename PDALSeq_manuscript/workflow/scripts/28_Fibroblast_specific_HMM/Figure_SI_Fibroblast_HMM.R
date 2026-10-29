library(tidyverse)
library(viridis)
library(cowplot)
library(ggpubr)

BASE <- "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10"
OUT  <- "results/28_Fibroblast_specific_HMM/Figure_SI_Fibroblast_HMM.svg"

Exp.levels <- c(
    "Ppaniscus_KB8711_PDALseq",  "Ppaniscus_KB14571_PDALseq",
    "Ggorilla_KB3781_PDALseq",   "Ggorilla_F7549_PDALseq",
    "Ppygmaeus_AG05252_PDALseq", "Pabelii_AG06213_PDALseq"
)
Exp.labels <- c(
    "P. paniscus KB8711",  "P. paniscus KB14571",
    "G. gorilla KB3781",   "G. gorilla F7549",
    "P. pygmaeus AG05252", "P. abelii AG06213"
)

NONB_levels <- c("APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS", "GQ")
NONB_labels <- c("APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS", "G4")
NONB.select <- c("APR", "CRU", "TRI", "Z", "DR", "G4")

ANN.select <- c("promoter_protCode_singleTrans", "utr5_protCode_singleTrans",
                "exon_protCode_singleTrans", "intron_protCode_singleTrans",
                "utr3_protCode_singleTrans", "enhancer_protCode_singleTrans",
                "repeats_protCode_singleTrans", "ngnr_protCode_singleTrans",
                "Recombination", "core_origins_hglft_genome_e0596_1cc9e0",
                "stochastic_origins_hglft_genome_2ff9ed_1ccfd0")
ANN.labels <- c("promoter", "5'-UTR", "exon", "intron", "3'-UTR",
                "enhancer", "repeats", "NGNR", "recombination",
                "core origins", "stochastic origins")

# ---- Order states by mean PDALseq signal (ascending = bottom of plot) ----
df.avg <- readr::read_csv(file.path(BASE, "Average_PDALseq_in_GHMM_states.csv"))
State.order <- df.avg %>%
    filter(Experiment %in% Exp.levels) %>%
    group_by(state) %>%
    summarize(grand_mean = mean(MEAN)) %>%
    arrange(grand_mean) %>%
    pull(state) %>%
    as.character()
# Label: rank-state_N (rank 1 = lowest signal). Update labels after biological interpretation.
State.labels <- paste0(seq_along(State.order), "-state_", State.order)

# ---- Panel A: mean PDALseq signal heatmap ----
df.avg$state      <- factor(as.character(df.avg$state), levels = State.order, labels = State.labels)
df.avg$Experiment <- factor(df.avg$Experiment, levels = Exp.levels, labels = Exp.labels)
df.avg            <- df.avg %>% filter(!is.na(Experiment))

P.tile <- ggplot(df.avg, aes(x = Experiment, y = state, fill = MEAN)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red", name = "Mean PDAL-seq signal") +
    theme_classic() +
    theme(axis.text.x   = element_text(angle = 25, hjust = 1, size = 6, color = "black"),
          axis.title.x  = element_blank(),
          legend.position  = "top",
          legend.title  = element_text(size = 8, color = "black"),
          legend.text   = element_text(size = 6, angle = 45, color = "black"),
          axis.title.y  = element_text(size = 6, color = "black"),
          legend.key.size = unit(0.1, "in"),
          axis.text.y   = element_text(size = 6, color = "black"),
          plot.margin   = margin(0, 0, 0, 0)) +
    ylab("Fibroblast HMM state")

# ---- Panel B: state size (Mbp) ----
df.count <- read.csv(file.path(BASE, "BP_in_state.csv"), header = FALSE)
colnames(df.count) <- c("state", "count")
df.count$state <- factor(as.character(df.count$state), levels = State.order, labels = State.labels)

P.count <- ggplot(df.count, aes(y = state, x = count / 1e6, fill = state,
                                label = round(count / 1e6))) +
    geom_bar(stat = "identity") +
    geom_text(size = 2) +
    theme_classic() +
    scale_fill_viridis_d(option = "plasma") +
    theme(axis.text.x   = element_text(angle = 90, size = 6, color = "black"),
          axis.title.x  = element_text(size = 8, color = "black"),
          legend.position  = "none",
          axis.title.y  = element_blank(),
          axis.text.y   = element_blank(),
          plot.title    = element_text(size = 8, color = "black"),
          plot.margin   = margin(0, 0, 0, 0)) +
    ggtitle("State size") +
    xlab("State size\n(Mbp)")

# ---- Panel C: conservation ----
df.con <- read.csv(file.path(BASE, "mean_state_conservation.csv"))
df.con$States <- factor(as.character(df.con$States), levels = State.order, labels = State.labels)

P.con <- ggplot(df.con, aes(y = States, x = "", fill = Mean_conservation)) +
    geom_tile() +
    theme_classic() +
    theme(axis.text.x   = element_blank(),
          axis.title.x  = element_text(size = 8, color = "black"),
          legend.position  = "top",
          legend.title  = element_text(size = 8, color = "black"),
          legend.text   = element_text(size = 6, color = "black", angle = 45, hjust = 1),
          axis.title.y  = element_blank(),
          axis.text.y   = element_blank(),
          legend.key.size = unit(0.1, "in"),
          plot.margin   = margin(0, 0, 0, 0)) +
    xlab("") +
    scale_fill_viridis("Conservation\nscore", trans = "log10")

# ---- Panel D: non-B enrichment ----
df.nonB <- read.csv(file.path(BASE, "State_nonB_enrichment", "compiled.csv"))
df.nonB$nonB  <- factor(df.nonB$nonB, levels = rev(NONB_levels), labels = rev(NONB_labels))
df.nonB$State <- factor(as.character(df.nonB$State), levels = State.order, labels = State.labels)
df.nonB$log2Enrichment[df.nonB$log2Enrichment >  log2(3)] <-  log2(3)
df.nonB$log2Enrichment[df.nonB$log2Enrichment < -log2(3)] <- -log2(3)
df.nonB$Sig   <- "ns"
df.nonB$Sig[df.nonB$Wilcox.p <= 0.05] <- "*"
df.nonB$Sig[df.nonB$Wilcox.p <= 0.01] <- "**"
df.nonB       <- df.nonB %>% filter(nonB %in% NONB.select)
df.nonB$LABEL <- paste(round(df.nonB$Enrichment, 1), df.nonB$Sig, sep = "\n")

P.nonB.E <- ggplot(df.nonB, aes(x = nonB, y = State, fill = log2Enrichment, label = LABEL)) +
    geom_tile() +
    geom_text(color = "black", size = 2) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limits = c(-log2(3), log2(3)),
                         name = "non-B DNA\nlog2Enrichment") +
    theme_classic() +
    theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.title.x  = element_blank(),
          legend.position  = "top",
          legend.title  = element_text(size = 8, color = "black"),
          legend.text   = element_text(size = 6, color = "black"),
          axis.title.y  = element_blank(),
          axis.text.y   = element_blank(),
          legend.key.size = unit(0.1, "in"),
          plot.margin   = margin(0, 0, 0, 0))

# ---- Panel E: functional annotation enrichment ----
df.ann <- read.csv(file.path(BASE, "State_functional_element_enrichment", "compiled.csv"))
df.ann$State <- factor(as.character(df.ann$State), levels = State.order, labels = State.labels)
df.ann$ann   <- factor(df.ann$ann, levels = ANN.select, labels = ANN.labels)
df.ann$log2Enrichment[df.ann$log2Enrichment >  log2(3)] <-  log2(3)
df.ann$log2Enrichment[df.ann$log2Enrichment < -log2(3)] <- -log2(3)
df.ann$Sig   <- "ns"
df.ann$Sig[df.ann$Wilcox.p <= 0.05] <- "*"
df.ann$Sig[df.ann$Wilcox.p <= 0.01] <- "**"
df.ann       <- df.ann %>% filter(!is.na(ann))
df.ann$LABEL <- paste(round(df.ann$Enrichment, 1), df.ann$Sig, sep = "\n")

P.ann.E <- ggplot(df.ann, aes(x = ann, y = State, fill = log2Enrichment, label = LABEL)) +
    geom_tile() +
    geom_text(color = "black", size = 2) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limits = c(-log2(3), log2(3)),
                         name = "Functional annotation\nlog2Enrichment") +
    theme_classic() +
    theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.title.x  = element_blank(),
          legend.position  = "top",
          legend.title  = element_text(size = 8, color = "black"),
          legend.text   = element_text(size = 6, color = "black"),
          axis.title.y  = element_blank(),
          axis.text.y   = element_blank(),
          legend.key.size = unit(0.1, "in"),
          plot.margin   = margin(0, 0, 0, 0))

# ---- Panel F: RNA expression ----
df.RNA <- read.csv(file.path(BASE, "State_RNA_expression", "compiled.csv"))
df.RNA$State <- factor(as.character(df.RNA$State), levels = State.order, labels = State.labels)

P.RNA <- ggplot(df.RNA, aes(y = State, x = Enrichment)) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = 1) +
    theme_classic() +
    theme(axis.text.x   = element_text(angle = 90, size = 6, color = "black"),
          axis.title.x  = element_text(size = 8, color = "black"),
          legend.position  = "none",
          axis.title.y  = element_blank(),
          axis.text.y   = element_blank(),
          plot.margin   = margin(0, 0, 0, 0)) +
    ggtitle("Transcription") +
    xlab("RNA-Seq\nenrichment")

# ---- Panel G: hypomethylated CpG ----
df.5mC  <- read.csv(file.path(BASE, "methylation", "compiled.csv"))
df.hypo <- df.5mC %>%
    filter(methylation_rate <= 0.2) %>%
    group_by(state) %>%
    summarize(Percent = sum(percent) * 100, .groups = "drop")
df.hypo$state <- factor(as.character(df.hypo$state), levels = State.order, labels = State.labels)

P.hypo <- ggplot(df.hypo, aes(x = Percent, y = state)) +
    geom_bar(stat = "identity", fill = "grey", color = "black") +
    theme_classic() +
    ggtitle("Hypo-methylated\nCpG sites") +
    xlab("Percent (%)") +
    theme(axis.text.x  = element_text(angle = 90, size = 6),
          axis.title.x = element_text(size = 8),
          legend.position = "none",
          axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          plot.margin  = margin(t = 0.23, r = 0, b = 0.175, l = 0, unit = "in"))

# ---- Assemble ----
P.A <- plot_grid(P.tile, P.con, P.count, P.nonB.E, P.ann.E, P.RNA, P.hypo,
                 align = "h", nrow = 1,
                 rel_widths = c(5, 1.2, 2, 5.2, 9.2, 2.5, 2.5),
                 labels = "A", label_size = 12)

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
ggsave(OUT, P.A, width = 6.5, height = 5, units = "in", scale = 1.0)
