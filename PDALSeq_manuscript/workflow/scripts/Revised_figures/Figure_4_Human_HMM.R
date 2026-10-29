
library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)

####A-B Correllation plots####
Key_size = 0.15

# NOTE: State.levels controls the y-axis order across all panels. After running
# step 27, examine data/27_Updated_Human_HMM_functional_analysis/Average_PDALseq_in_GHMM_states.csv
# and re-order states from lowest to highest mean PDALseq signal across experiments.
State.levels = c("2", "4", "6", "1", "0", "5", "8", "7", "9", "3")
State.labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

Exp.levels = c("Hsapien_GM24385_PDALseq",
               "Hsapien_GM24143_PDALseq",
               "Hsapien_HEK293T_PDALseq",
               "Hsapien_CCL86-Raji_PDALseq",
               "Hsapien_K562_PDALseq")

Exp.label = c("GM24385-HG002",
              "GM24143-HG004",
              "HEK293T",
              "CCL86-Raji",
              "K562")

NONB_levels = c("all", "APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS", "GQ")
NONB_labels = c("all", "APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS", "G4")

####Hidden Markov model#####

df = readr::read_csv("data/27_Updated_Human_HMM_functional_analysis/Average_PDALseq_in_GHMM_states.csv") 

df %>% group_by(state) %>% summarize(MEAN = mean(MEAN)) %>% arrange(MEAN)

df$state      = factor(df$state, levels = State.levels)
df$Experiment = factor(df$Experiment, levels = Exp.levels, labels = Exp.label)

P.tile = ggplot(df, aes(x = Experiment, y = state, fill = MEAN)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red", name = "Mean PDAL-seq signal") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.title.x = element_blank(),
          legend.position = "top",
          legend.title = element_text(size = 8, color = "black"),
          legend.text = element_text(size = 6, color = "black"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          legend.key.size = unit(0.1, "in"))


####Make a plot that shows the number of blocks in each state####

Human_HMM_color = c("#386350", "#008A69", "#00C992", "#B6DBFF", "#7BB0DF", "#1964B0", "#E9DC6D", "#894B45", "#F4A637", "#D11F5E")

df.count = read.csv("data/27_Updated_Human_HMM_functional_analysis/BP_in_state.csv", header = F)
colnames(df.count) = c("state", "count")
df.count$state   = factor(df.count$state, levels = State.levels, labels = State.labels)

P.count = ggplot(df.count, aes(y = state, x = count/1000000, fill = state, label = round(count/1000000))) +
    geom_bar(stat = "identity") +
    geom_text(size = 2) +
    theme_classic() +
    scale_fill_manual(values = Human_HMM_color) +
    theme(axis.text.x = element_text(angle = 90, size = 6, color = "black"),
          axis.title.x = element_text(size = 8, color = "black"),
          legend.position = "none",
          axis.title.y = element_text(size = 8, color = "black"),
          axis.text.y  = element_text(size = 6, color = "black"),
          plot.title   = element_text(size = 8, color = "black"),
          plot.margin  = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          legend.key.size = unit(0.1, "in")) +
    xlab("Size (Mbp)") +
    ylab("Humman HMM State")

####make a plot that shows non-B motif enrichment in each state####

df.nonB = read.csv("data/27_Updated_Human_HMM_functional_analysis/State_nonB_enrichment/hsapiens.csv")

NONB.select = c("APR", "CRU", "TRI", "Z", "DR", "G4")
df.nonB$log2Enrichment[which(df.nonB$log2Enrichment >  log2(3))] =  log2(3)
df.nonB$log2Enrichment[which(df.nonB$log2Enrichment < -log2(3))] = -log2(3)
df.nonB$nonB  = factor(df.nonB$nonB,  levels = rev(NONB_levels), labels = rev(NONB_labels))
df.nonB$State = factor(df.nonB$State, levels = State.levels)
df.nonB$Sig   = rep("ns", nrow(df.nonB))
df.nonB$Sig[which(df.nonB$Wilcox.p <= 0.05)] = "*"
df.nonB$Sig[which(df.nonB$Wilcox.p <= 0.01)] = "**"
df.nonB = df.nonB %>% filter(nonB %in% NONB.select)
df.nonB$LABEL = paste(round(df.nonB$Enrichment, digits = 1), df.nonB$Sig, sep = "\n")

P.nonB.E = ggplot(df.nonB, aes(x = nonB, y = State, fill = log2Enrichment, label = LABEL)) +
    geom_tile() +
    geom_text(color = "black", size = 2) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-log2(3), log2(3)),
                         name = "non-B DNA motif\nlog2Enrichment") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.title.x = element_blank(),
          legend.position = "top",
          legend.title = element_text(size = 8, color = "black"),
          legend.text  = element_text(size = 6, color = "black"),
          axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          legend.key.size = unit(0.1, "in"),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

####make a plot that shows functional motif enrichment in each state####

df.ann = read.csv("data/27_Updated_Human_HMM_functional_analysis/State_functional_element_enrichment/hsapiens.csv")

ANNOTATIONS = c("promoter_protCode_singleTrans", "utr5_protCode_singleTrans",
                "gene_protCode_singleTrans", "gene_nonprotCode",
                "cds_protCode_singleTrans", "exon_protCode_singleTrans",
                "intron_protCode_singleTrans", "utr3_protCode_singleTrans",
                "enhancer_protCode_singleTrans", "rna_protCode_singleTrans",
                "chrG_cpgislands", "core_stochastic_origins_chm13v2.0_hglft.ucsc",
                "core_origins_hglft_genome_e0596_1cc9e0",
                "stochastic_origins_hglft_genome_2ff9ed_1ccfd0",
                "Recombination", "repeats_protCode_singleTrans",
                "ngnr_protCode_singleTrans")

ANNOTATIONS.labels = c("promoter", "5'-UTR", "gene", "gene not protein coding",
                       "cds", "exon", "intron", "3'-UTR", "enhancer", "RNA",
                       "CpG islands", "replication origins", "core origins",
                       "stochastic origins", "recombination sites", "repeats", "NGNR")

df.ann$State = factor(df.ann$State, levels = State.levels)
df.ann$ann   = factor(df.ann$ann,   levels = ANNOTATIONS, labels = ANNOTATIONS.labels)
df.ann$log2Enrichment[which(df.ann$log2Enrichment >  log2(3))] =  log2(3)
df.ann$log2Enrichment[which(df.ann$log2Enrichment < -log2(3))] = -log2(3)
df.ann$Sig = rep("ns", nrow(df.ann))
df.ann$Sig[which(df.ann$Wilcox.p <= 0.05)] = "*"
df.ann$Sig[which(df.ann$Wilcox.p <= 0.01)] = "**"
df.ann = df.ann %>%
    filter(ann %in% c("promoter", "enhancer", "5'-UTR", "exon", "intron", "3'-UTR",
                      "repeats", "NGNR", "recombination sites",
                      "core origins", "stochastic origins"))
df.ann$LABEL = paste(round(df.ann$Enrichment, digits = 1), df.ann$Sig, sep = "\n")

P.ann.E = ggplot(df.ann, aes(x = ann, y = State, fill = log2Enrichment, label = LABEL)) +
    geom_tile() +
    geom_text(color = "black", size = 2) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-log2(3), log2(3)),
                         name = "Functional annotation\nlog2Enrichment") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.title.x = element_blank(),
          legend.position = "top",
          legend.title = element_text(size = 8, color = "black"),
          legend.text  = element_text(size = 6, color = "black"),
          axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          legend.key.size = unit(0.1, "in"),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))


####Make a plot that shows human HMM states intersecting with ape HMM states####

df.intersect = read.csv("data/27_Updated_Human_HMM_functional_analysis/intersect_results.csv")
l.human_state = strsplit(df.intersect$query,    split = "/")
l.ape_state   = strsplit(df.intersect$database, split = "/")
human_state = c()
ape_state   = c()

for (i in 1:nrow(df.intersect)) {
    human_state[i] = gsub(".bed", "", gsub("state_", "", l.human_state[[i]][3]))
    ape_state[i]   = gsub(".bed", "", gsub("state_", "", l.ape_state[[i]][6]))
}

df.intersect$human_state = factor(human_state, levels = State.levels)
df.intersect$ape_state   = factor(ape_state,   levels =  c("2", "5", "0", "6", "1", "3", "7", "4"), labels = rev(c("8-origin", "7-recomb", "6-gene", "5-variable", "4-variable", "3-low", "2-low", "1-low")))
df.intersect$Percent     = 100 * df.intersect$intersect / df.intersect$total

HUMAN_states = unique(df.intersect$human_state)
Percent = c()

for (i in 1:length(HUMAN_states)) {
    df.i       = df.intersect %>% filter(human_state == HUMAN_states[i])
    Percent[i] = 100 - sum(df.i$Percent)
}

df.intersect1 = bind_rows(
    df.intersect %>% select(human_state, ape_state, Percent),
    data.frame("human_state" = factor(HUMAN_states, levels = State.levels),
               ape_state = NA, "Percent" = Percent)
)

P.intersect = ggplot(df.intersect1, aes(human_state, x = Percent, fill = ape_state)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    scale_fill_manual(name = "Ape state\nintersect",
                      values = c("#386350", "#386350", "#386350", "#B6DBFF", "#B6DBFF",
                                 "#882D71", "#DB5829", "#F4A637")) +
    xlab("Percent human state\nin ape states (%)") +
    scale_x_continuous(limits = c(0, 100)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.title.x = element_text(size = 8, color = "black"),
          legend.position = "top",
          legend.title = element_text(size = 8, color = "black"),
          legend.text  = element_text(size = 6, color = "black"),
          axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          legend.key.size = unit(0.1, "in"),
          plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "in"))

####Human states mapped to chromosomes####

df.bed   = readr::read_delim(
    "data/09_Human_HMM/01_Remove_sex_chromosomes/Hsapien_GM24385_PDALseq.bg",
    delim = "\t", col_names = c("chr", "start", "end", "reads"))[, 1:3]
df.state = readr::read_csv(
    "data/26_Ape_PCA_and_Human_HMM_with_K562/Human_GHMM/10/GHMM_states.csv")
df.bed$hidden_state = df.state$hidden_state
df.bed$hidden_state = factor(df.bed$hidden_state, levels = State.levels, labels = State.labels)

CHR = unique(df.bed$chr)
df.bed$chr = factor(df.bed$chr, levels = CHR)

df.chr.sizes = readr::read_delim("resources/genomes/Hsapien.chrom.sizes",
                                  col_names = c("chr", "end"))
df.chr.sizes$chr = factor(df.chr.sizes$chr, levels = CHR)
df.chr.sizes = df.chr.sizes %>% filter(!is.na(chr))

df.cen = readr::read_delim("resources/centromeres/Hsapien.bed",
                             col_names = c("chr", "start", "end", "ann"))
df.cen$chr = factor(df.cen$chr, levels = CHR)
df.cen = df.cen %>% filter(!is.na(chr))

l.df.bed = df.bed %>% group_by(chr) %>% group_split()

l.df = list()

for (i in 1:length(l.df.bed)) {
    df.i = data.frame(l.df.bed[[i]] %>%
                          group_by(hidden_state) %>%
                          summarize(BP = sum(end - start + 1)))
    Size = df.chr.sizes$end[which(df.chr.sizes$chr == l.df.bed[[i]]$chr[1])]
    df.i.nonsyntenic = data.frame("hidden_state" = "blacklist",
                                  "BP" = Size - sum(df.i$BP))
    df.i = bind_rows(df.i, df.i.nonsyntenic)
    df.i$Percent = 100 * df.i$BP / Size
    df.i$chr = l.df.bed[[i]]$chr[1]
    l.df[[i]] = df.i
}

df.state.percent = bind_rows(l.df)
df.state.percent$chr = factor(df.state.percent$chr, levels = CHR)
df.state.percent$hidden_state = factor(df.state.percent$hidden_state, levels = c("blacklist", State.levels))
P.ideogram = ggplot() +
    geom_rect(data = df.bed,
              mapping = aes(xmin = -0.8, xmax = 0.2, ymin = start, ymax = end,
                            fill = hidden_state)) +
    scale_fill_manual(values = c("grey", Human_HMM_color)) +
    geom_rect(data = df.chr.sizes,
              mapping = aes(xmin = -1.5, xmax = -1.0, ymin = 0, ymax = end),
              fill = "grey", color = "black") +
    geom_rect(data = df.cen,
              mapping = aes(xmin = -1.6, xmax = -0.9, ymin = start, ymax = end),
              fill = "dimgrey", color = "black") +
    facet_wrap(~chr, nrow = 1) +
    theme_void() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          strip.text.y = element_blank(),
          strip.background = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x  = element_blank(),
          panel.spacing = unit(0, 'lines'))

P.percent = ggplot(df.state.percent, aes(x = chr, y = Percent, fill = hidden_state)) +
    facet_wrap(~chr, nrow = 1, scales = "free_x") +
    ggtitle("Human state by chromosome") +
    geom_bar(stat = "sum") +
    theme_classic() +
    scale_fill_manual(values = c("grey", Human_HMM_color)) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text.x  = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "top",
          axis.title.y = element_text(size = 8, color = "black"),
          axis.text.y  = element_text(size = 6, color = "black"),
          plot.title   = element_text(size = 8, color = "black"),
          strip.text   = element_text(size = 6, color = "black"),
          panel.spacing = unit(0, 'lines')) +
    ylab("Percent (%)")

####RM family enrichment in states####

df.states    = read.csv("data/27_Updated_Human_HMM_functional_analysis/Human_state_rm_intersect_bp.csv")
df.baseline  = read.csv("data/12_PDALseq/baseline_human_state_rm_intersect_bp.csv")
df.intersect = bind_rows(df.states, df.baseline)
df.intersect$Percent = 100 * df.intersect$Intersect_bp / df.intersect$Total_bp_state

l.Class = strsplit(as.character(df.intersect$Ann), split = "/")

Parent_class = c()
for (i in 1:nrow(df.intersect)) {
    Parent_class[i] = l.Class[[i]][1]
}

Parent_class_levels = c("Unknown", "Low_complexity",
                         "Satellite", "Simple_repeat",
                         "tRNA", "srpRNA", "snRNA", "scRNA", "rRNA",
                         "DNA?", "DNA", "RC",
                         "LTR", "Retroposon", "SINE", "LINE")
Parent_class_labels = c("Unknown", "Low_complexity",
                         "Satellite", "Simple_repeat",
                         "tRNA", "srpRNA", "snRNA", "scRNA", "rRNA",
                         "DNA?", "DNA transposon", "RC",
                         "LTR", "Retroposon", "SINE", "LINE")

df.intersect$Parent_class = factor(Parent_class,
                                    levels = Parent_class_levels,
                                    labels = Parent_class_labels)

df.intersect$State = factor(
    as.character(df.intersect$State),
    levels = c("Human HMM Region", paste("state_", State.levels, sep = "")),
    labels = c("Human HMM\nRegion", State.labels)
)

P.RM.intersect = ggplot(df.intersect, aes(y = State, x = Percent, fill = State)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    geom_vline(data = df.intersect %>%
                   filter(State == "Human HMM\nRegion") %>%
                   summarize(SUM = sum(Percent)),
               mapping = aes(xintercept = SUM)) +
    scale_fill_manual(values = c("black", Human_HMM_color)) +
    ggtitle("Repeat content of\nhuman states") +
    theme(axis.text.x  = element_text(size = 6, color = "black"),
          axis.title.x = element_text(size = 8, color = "black"),
          legend.position = "none",
          axis.title.y = element_text(size = 8, color = "black"),
          axis.text.y  = element_text(size = 6, color = "black"),
          plot.title   = element_text(size = 8, color = "black"),
          strip.text   = element_text(size = 8, color = "black"),
          plot.margin  = margin(t = 0, unit = "in")) +
    xlab("Percent (%)") +
    ylab("Hidden Markov State")

Exclude = c("Unknown", "Low_complexity", "tRNA", "srpRNA", "snRNA", "scRNA", "rRNA",
            "DNA?", "RC", "Retroposon")

df.parC = df.intersect %>%
    filter(!Parent_class %in% Exclude) %>%
    arrange(Parent_class)
df.parC$Parent_class = factor(df.parC$Parent_class,
                               levels = unique(df.parC$Parent_class),
                               labels = c("Satellite", "Simple repeat",
                                          "DNA transposon", "LTR", "SINE", "LINE"))

df.rm_E = read.csv("results/27_Updated_Human_HMM_functional_analysis/RM_family_enrichment.csv")

df.rm_E$State = factor(df.rm_E$State,
    levels = State.levels,
    labels =  State.labels) 
df.rm_E$Sig = rep("ns", nrow(df.rm_E))
df.rm_E$Sig[which(df.rm_E$Wilcox.p <= 0.05)] = "*"
df.rm_E$Sig[which(df.rm_E$Wilcox.p <= 0.01)] = "**"
df.rm_E$LABEL = paste(round(df.rm_E$Enrichment, digits = 1), df.rm_E$Sig, sep = " ")
df.rm_E = df.rm_E %>% filter(!Ann %in% c("Psuedogene", "Low_complexity"))
df.rm_E$Ann = factor(df.rm_E$Ann,
                      levels = c("Satelitte", "Simple_repeat", "DNA_transposons",
                                 "LTR", "SINE", "LINE"),
                      labels = c("Satellite", "Simple repeat", "DNA transposon",
                                 "LTR", "SINE", "LINE"))

l.df = df.parC %>% group_by(State) %>% group_split()
l.df.out = list()
for (i in 1:length(l.df)) {
    df.i = l.df[[i]] %>%
        group_by(Parent_class) %>%
        summarise(Percent = sum(Percent, na.rm = TRUE))
    df.i$State  = l.df[[i]]$State[1]
    l.df.out[[i]] = df.i
}

df.out = bind_rows(l.df.out)
LABEL  = c()
df.rm_E$Parent_class = df.rm_E$Ann

for (i in 1:nrow(df.out)) {
    df.i = df.rm_E %>%
        filter(Ann == df.out$Parent_class[i]) %>%
        filter(State == as.character(df.out$State[i]))
    LABEL[i] = if (nrow(df.i) == 0) NA else df.i$LABEL[1]
}
df.out$LABEL = LABEL

df.all = df.intersect
df.all$Parent_class = "all"
df.all = bind_rows(df.all, df.parC %>% filter(Parent_class != "Satellite"))
df.all$Parent_class = factor(
    df.all$Parent_class,
    levels = c("all", "Simple repeat", "DNA transposon", "LTR", "SINE", "LINE"),
    labels = c("all", "simple repeat", "DNA transposon", "LTR", "SINE", "LINE")
)

df.out = df.out %>% filter(Parent_class != "Satellite")
df.out$Parent_class = factor(
    df.out$Parent_class,
    levels = c("all", "Simple repeat", "DNA transposon", "LTR", "SINE", "LINE"),
    labels = c("all", "simple repeat", "DNA transposon", "LTR", "SINE", "LINE")
)

P.RM.class = ggplot(df.all, aes(y = State, x = Percent, fill = State)) +
    geom_bar(stat = "identity") +
    facet_wrap(~Parent_class, nrow = 1, scales = "free_x") +
    geom_text(data = df.out, mapping = aes(x = 30, label = LABEL), size = 2) +
    theme_classic() +
    scale_fill_manual(values = c("black", Human_HMM_color)) +
    geom_vline(data = df.all %>%
                   filter(State == "Human HMM\nRegion") %>%
                   group_by(Parent_class) %>%
                   summarize(SUM = sum(Percent)),
               mapping = aes(xintercept = SUM),
               color = "dimgrey", linetype = 'dotted', size = 0.5) +
    xlab("Percent of human HMM state in repeats (%)") +
    theme(axis.text.x  = element_text(size = 6, color = "black"),
          axis.title.x = element_text(size = 8, color = "black"),
          legend.position = "none",
          axis.title.y = element_text(size = 8, color = "black"),
          axis.text.y  = element_text(size = 6, color = "black"),
          plot.title   = element_text(size = 8, color = "black"),
          strip.text   = element_text(size = 6, color = "black"),
          plot.margin  = margin(t = 0.2, r = 0.5, b = 0, l = 0.5, unit = "in")) +
  ylab("Human HMM state")

####Compile plots####

P.A = plot_grid(P.count, P.tile, P.nonB.E, P.ann.E, P.intersect,
                align = "h", nrow = 1,
                rel_widths = c(3, 5, 5.1, 7.9, 6.8),
                labels = c("A"), label_size = 12)

plot_grid(P.A, P.RM.class, ncol = 1)

P = plot_grid(P.A, P.percent, P.percent, P.RM.class,
              ncol = 1, rel_heights = c(1.0, 0.4, 0.3, 0.5),
              labels = c("", "B", "C"), label_size = 12)

ggsave("revised_plots/Figure_4.svg",          P,           width = 6.5, height = 7,  units = "in")
ggsave("revised_plots/Figure_4_ideogram.png", P.ideogram,  width = 6.5, height = 2,  units = "in", dpi = 1200)
