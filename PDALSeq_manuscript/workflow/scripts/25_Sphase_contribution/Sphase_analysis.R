library(tidyverse)
library(cowplot)
library(matrixTests)
library(viridis)
library(ggpubr)

# ── Sample metadata ───────────────────────────────────────────────────────────

dmso_GM24143 <- c(
    "Hsa_GM24143_DMSO_PDALSeq_js4012d20mM_S4_R1_001",
    "Hsa_GM24143_DMSO_PDALSeq_js4012e20mM_S5_R1_001",
    "Hsa_GM24143_DMSO_PDALSeq_js4012f20mM_S6_R1_001"
)
dmso_GM24385 <- c(
    "Hsa_GM24385_DMSO_PDALSeq_js4012a20mM_S1_R1_001",
    "Hsa_GM24385_DMSO_PDALSeq_js4012b20mM_S2_R1_001",
    "Hsa_GM24385_DMSO_PDALSeq_js4012c20mM_S3_R1_001"
)
bulk_GM24143 <- c(
    "Hsapien_GM24143_LCL_PDALseq_js2007d20mM_S2_R1_001",
    "Hsapien_GM24143_LCL_PDALseq_js2007e20mM_S3_R1_001",
    "Hsapien_GM24143_LCL_PDALseq_js2007f40mM_S4_R1_001",
    "Hsapien_GM24143_LCL_PDALseq_js2007g40mM_S5_R1_001"
)
bulk_GM24385 <- c(
    "Hsapien_GM24385_LCL_PDALseq_js1014A20_1_S2_R1_001",
    "Hsapien_GM24385_LCL_PDALseq_js1014A20_2_S3_R1_001",
    "Hsapien_GM24385_LCL_PDALseq_js1014A40_1_S4_R1_001",
    "Hsapien_GM24385_LCL_PDALseq_js1014A40_2_S5_R1_001"
)

all_samples <- c(dmso_GM24143, dmso_GM24385, bulk_GM24143, bulk_GM24385)

# ── Load data ─────────────────────────────────────────────────────────────────

df <- readr::read_csv("results/25_Sphase_contribution/Sphase_1kb_coverage.csv")

nrow(df)

# ── Section 1: QC ─────────────────────────────────────────────────────────────

# Remove sex chromosomes and unplaced contigs
df <- df %>% filter(!chr %in% c("chrX", "chrY", "chrM"),
                    grepl("^chr[0-9]", chr))
nrow(df)

# Summary stats per sample
for (s in all_samples) {
    x <- df[[s]]
    cat(sprintf("%-65s  min=%6.2f  median=%6.2f  mean=%6.2f  max=%8.2f\n",
                s, min(x, na.rm=TRUE), median(x, na.rm=TRUE),
                mean(x, na.rm=TRUE), max(x, na.rm=TRUE)))
}

# Histograms of RPM distribution by condition
df_long <- df %>%
    pivot_longer(cols = all_of(all_samples),
                 names_to = "Sample", values_to = "RPM") %>%
    mutate(Condition = case_when(
        Sample %in% c(dmso_GM24143, dmso_GM24385) ~ "DMSO",
        TRUE ~ "Bulk"
    ),
    Cell_line = case_when(
        Sample %in% c(dmso_GM24143, bulk_GM24143) ~ "GM24143",
        TRUE ~ "GM24385"
    ))

P.hist <- ggplot(df_long, aes(x = log1p(RPM), fill = Condition)) +
    geom_histogram(bins = 100, alpha = 0.6, position = "identity") +
    facet_wrap(~ Cell_line) +
    theme_classic() +
    xlab("log1p(RPM)") +
    scale_fill_manual(values = c("Bulk" = "#2196F3", "DMSO" = "#F44336"))

# Within-condition correlation matrix (spot-check tech-rep concordance)
cor_dmso  <- cor(df[, c(dmso_GM24143, dmso_GM24385)],  method = "pearson")
cor_bulk  <- cor(df[, c(bulk_GM24143,  bulk_GM24385)],  method = "pearson")
cat("\nDMSO Pearson correlations:\n"); print(round(cor_dmso, 3))
cat("\nBulk Pearson correlations:\n"); print(round(cor_bulk, 3))

# ── Section 2: log2FC DMSO vs Bulk per cell line ──────────────────────────────

pseudocount <- 0.001   # on RPM scale; adjust if needed

# Vectorized row standard error — avoids slow apply() over 3M rows
row_se <- function(mat) {
    n     <- ncol(mat)
    means <- rowMeans(mat, na.rm = TRUE)
    sd_v  <- sqrt(rowSums((mat - means)^2, na.rm = TRUE) / (n - 1))
    sd_v / sqrt(n)
}

compute_fc <- function(df, dmso_cols, bulk_cols, label) {
    dmso_mat  <- as.matrix(df[, dmso_cols])
    bulk_mat  <- as.matrix(df[, bulk_cols])
    mean_dmso <- rowMeans(dmso_mat, na.rm = TRUE)
    mean_bulk <- rowMeans(bulk_mat, na.rm = TRUE)
    se_dmso   <- row_se(dmso_mat)
    se_bulk   <- row_se(bulk_mat)
    log2fc    <- log2((mean_bulk + pseudocount) / (mean_dmso + pseudocount))
    se_log2fc <- (1 / log(2)) * sqrt((se_bulk  / (mean_bulk  + pseudocount))^2 +
                                     (se_dmso  / (mean_dmso  + pseudocount))^2)
    z         <- (log2fc - mean(log2fc, na.rm = TRUE)) / sd(log2fc, na.rm = TRUE)
    wt        <- row_wilcoxon_twosample(dmso_mat, bulk_mat)
    tt        <- row_t_welch(dmso_mat, bulk_mat)
    data.frame(chr        = df$chr,
               start      = df$start,
               end        = df$end,
               mean_dmso  = mean_dmso,
               se_dmso    = se_dmso,
               mean_bulk  = mean_bulk,
               se_bulk    = se_bulk,
               log2FC     = log2fc,
               se_log2fc  = se_log2fc,
               z_score    = z,
               wilcox_p   = wt$pvalue,
               t_p        = tt$pvalue,
               Cell_line  = label)
}

fc_143 <- compute_fc(df, dmso_GM24143, bulk_GM24143, "GM24143")
fc_385 <- compute_fc(df, dmso_GM24385, bulk_GM24385, "GM24385")

# Distribution of log2FC

df_fc = bind_rows(fc_143, fc_385)

P.fc <- df_fc %>%
    ggplot(aes(x = log2FC, colour = Cell_line)) +
    geom_density() +
    theme_classic() +
    geom_vline(xintercept = 0, linetype = "dashed")

P.Zscore <- df_fc %>%
  ggplot(aes(x =  z_score, colour = Cell_line)) +
  geom_density() +
  theme_classic() +
  geom_vline(xintercept = c(-2,2)) +
  geom_vline(xintercept = 0, linetype = "dashed")

length(which(df_fc$wilcox_p <= 0.05))
df_fc$neglog10wilcox_P = -log10(df_fc$wilcox_p)

df_comapare = bind_cols(df_fc %>% filter(Cell_line == "GM24143") %>% select(chr, start, end, log2FC), df_fc %>% filter(Cell_line != "GM24143") %>% select(log2FC))

P.Cell.Cor = df_comapare %>% ggplot(aes(x =  log2FC...4, y = log2FC...5)) +
  geom_hex() +
  geom_abline(slope = 1, intercept = 1, color = "red") +
  theme_classic() +
  stat_cor(method = "pearson") +
  xlab("log2FC GM24143") +
  ylab("log2FC GM24385") +
  scale_fill_viridis(trans = "log10") 

P.volcano = df_fc %>%
  ggplot(aes(x =  log2FC, y = t_p)) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 0.05, color = "red") +
  theme_classic() +
  facet_wrap(~Cell_line) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_continuous(scales::transform_compose("log10", "reverse"))

df_count = data.frame(table(df_fc$t_p))
head(df_count %>% arrange(-Freq))

df_fc$Sig = "ns"
df_fc$Sig[which(df_fc$t_p <= 0.05)] = "p < 0.05"
df_fc$Sig[which(df_fc$t_p <= 0.01)] = "p < 0.01"
df_fc$Sig[which(df_fc$t_p <= 0.001)] = "p < 0.001"
df_fc$Sig = factor(df_fc$Sig, levels = c("ns", "p < 0.05", "p < 0.01", "p < 0.001"))

df_fc %>% filter(chr == "chr9") %>% ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = log2FC, fill = Sig)) +
  geom_rect(stat = "identity") +
  facet_wrap(~Cell_line, ncol = 1) + 
  theme_classic() +
  scale_fill_manual(values = c("grey", "black", "orange", "red", "darkred")) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "red")



df_fc %>% filter(chr == "chr1") %>% ggplot(aes(x = start, y =log2FC, color = Sig)) +
  geom_point(size = 0.1) +
  facet_wrap(~Cell_line, ncol = 1) + 
  theme_classic() +
  scale_color_manual(values = c("grey", "black", "orange", "red", "darkred")) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "red") +
  geom_hline(yintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, color = "black") 
# ── Section 3: Identify S-phase-associated regions ───────────────────────────

threshold <- log2(1.5)   # Fold change cutoff; adjust after reviewing the distribution

head(df_fc)
df_fc$PDALSeq_increase = "none"

# Regions with higher signal bulk phase (S-phase signal)

df_fc$PDALSeq_increase[which(df_fc$log2FC >= threshold)] = "S-phase"

# Regions with Higher signal in DMSO (G1/G0 arrest)

df_fc$PDALSeq_increase[which(df_fc$log2FC <= -threshold)] = "G0/G1 arrest"

df_fc$sum = df_fc$end - df_fc$start + 1

P_sum = df_fc %>% ggplot(aes(x = Cell_line, y = sum, fill = PDALSeq_increase)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = viridis(3))

df_chr9 = df_fc %>% filter(chr == "chr9") %>% filter(Cell_line == "GM24385")
XLIM = range(df_chr9$start)

P_9_FC =  df_chr9 %>% ggplot(aes(x = start, y =log2FC, color = Sig)) +
  geom_point(size = 0.1) +
  theme_classic() +
  scale_color_manual(values = c("grey", "black", "orange", "red", "darkred")) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "red") +
  geom_hline(yintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, color = "black") +
  scale_x_continuous(limits = XLIM, expand = c(0,0)) 

P_9_FC

P_9_bin =  df_chr9 %>% ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = PDALSeq_increase)) +
  geom_rect() +
  theme_void() +
  scale_fill_manual(values = viridis(3)) +
  scale_x_continuous(limits = XLIM, expand = c(0,0)) 

P_compare = plot_grid(P_9_FC, P_9_bin, rel_heights = c(3,1), align = "v", ncol = 1)


# ── Section 4: Write BED files and processed data for plotting─────────────────

# Sort before writing; consider running `bedtools merge -i` on the output
write_bed <- function(df, path) {
    df %>%
        arrange(chr, start) %>%
        select(chr, start, end) %>%
        mutate(start = as.integer(start), end = as.integer(end)) %>%
        write.table(path, sep = "\t", quote = FALSE,
                    row.names = FALSE, col.names = FALSE)
    cat("Written:", path, "\n")
}

unique(df_fc$PDALSeq_increase)
unique(df_fc$Cell_line)
write_bed(df_fc %>% filter(PDALSeq_increase == "S-phase") %>% filter(Cell_line == "GM24385"),
          "results/25_Sphase_contribution/GM24385/Sphase_increases_ssDNA_windows.bed")
write_bed(df_fc %>% filter(PDALSeq_increase == "G0/G1 arrest") %>% filter(Cell_line == "GM24385"),
          "results/25_Sphase_contribution/GM24385/G0_G1_arrest_increases_ssDNA_windows.bed")
write_bed(df_fc %>% filter(PDALSeq_increase == "S-phase") %>% filter(Cell_line == "GM24143"),
          "results/25_Sphase_contribution/GM24143/Sphase_increases_ssDNA_windows.bed")
write_bed(df_fc %>% filter(PDALSeq_increase == "G0/G1 arrest") %>% filter(Cell_line == "GM24143"),
          "results/25_Sphase_contribution/GM24143/G0_G1_arrest_increases_ssDNA_windows.bed")

# Merge adjacent windows (run in shell after reviewing):
# bedtools merge -i results/25_Sphase_contribution/Sphase_enriched_windows.bed \
#   > results/25_Sphase_contribution/Sphase_enriched_merged.bed

# Write the data to a new file for plotting

write.csv(df_fc, "results/25_Sphase_contribution/Cell_cycle_effects_on_PDALSeq_1kbp_windows.csv", row.names = F)
