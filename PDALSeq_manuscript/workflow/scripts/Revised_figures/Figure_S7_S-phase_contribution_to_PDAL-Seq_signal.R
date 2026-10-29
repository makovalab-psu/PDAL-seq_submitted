library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)

####Theme set####

theme_set(theme_classic(base_size = 6) +
            theme(text = element_text(color = "black", size = 6),
                  axis.text = element_text(color = "black", size = 6),
                  axis.title = element_text(color = "black", size = 6),
                  strip.text = element_text(color = "black", size = 6),
                  legend.text = element_text(color = "black", size = 6),
                  legend.title = element_text(color = "black", size = 6),
                  plot.title = element_text(color = "black", size = 6)))

####Correlation####

df = readr::read_csv("results/25_Sphase_contribution/Technical_replicate_correlation.csv") %>% filter(Window_size == 1000000)

Sample_levels = unique(df$Sample_1)
df$Sample_1 = factor(df$Sample_1, levels = Sample_levels)
df$Sample_2 = factor(df$Sample_2, levels = Sample_levels)

P = ggplot(df, aes(x = Sample_1, y = Sample_2, fill = Spearman)) +
	geom_tile() +
	theme_classic() +
	theme(axis.text.x = element_text(angle =90, size = 2),
	      axis.text.y = element_blank(),
	      axis.title = element_blank()) +
	scale_fill_viridis(limits = c(0,1), option = "turbo") +
	coord_fixed()

####FC will cell cycle arrest results####

df_fc = readr::read_csv("results/25_Sphase_contribution/Cell_cycle_effects_on_PDALSeq_1kbp_windows.csv")

unique(df_fc$Cell_line)

df_chr9 = df_fc %>% filter(chr == "chr9")

XLIM = range(df_chr9$start)
XLIM[1] = XLIM[1] - 1
XLIM[2] = XLIM[2] + 1

P_9_FC_GM24385 =  df_chr9  %>% filter(Cell_line == "GM24385") %>% ggplot(aes(x = start, y =log2FC, color = Sig)) +
  geom_point(size = 0.1) +
  theme_classic() +
  scale_color_manual(values = c("grey", "black", "orange", "red", "darkred")) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "red") +
  geom_hline(yintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, color = "black") +
  scale_x_continuous(limits = XLIM, expand = c(0,0)) 

P_9_FC_GM24143 =  df_chr9  %>% filter(Cell_line == "GM24143") %>% ggplot(aes(x = start, y =log2FC, color = Sig)) +
  geom_point(size = 0.1) +
  theme_classic() +
  scale_color_manual(values = c("grey", "black", "orange", "red", "darkred")) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "red") +
  geom_hline(yintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, color = "black") +
  scale_x_continuous(limits = XLIM, expand = c(0,0)) 
df_chr9$PDALSeq_increase  = factor(df_chr9$PDALSeq_increase, levels = c("S-phase", "G0/G1 arrest", "none"), labels = c("S-phase", "G0/G1 arrest", "neither"))

P_9_bin =  ggplot()  + 
  geom_rect(data = df_chr9  %>% filter(Cell_line == "GM24385") %>% filter(PDALSeq_increase != "none"),
            aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = PDALSeq_increase)) +
  geom_rect(data = df_chr9  %>% filter(Cell_line == "GM24143") %>% filter(PDALSeq_increase != "none"),
            aes(xmin = start, xmax = end, ymin = 1, ymax = 2, fill = PDALSeq_increase)) +
  theme_void() +
  scale_fill_manual(values = c("#E95562FF",  "#331068FF", "grey")) +
  scale_x_continuous(limits = XLIM, expand = c(0,0)) 	

P_FC = plot_grid(P_9_FC_GM24385, P_9_FC_GM24143, P_9_bin, rel_heights = c(1, 1, 0.1), ncol = 1, align = "v")

####Fraction in each cell cycle association####

df_sum_GM24385 = df_fc %>% filter(Cell_line == "GM24385") %>% group_by(PDALSeq_increase) %>% summarise(SUM = sum(sum))
df_sum_GM24385$Cell_line = "GM24385"
df_sum_GM24143 = df_fc %>% filter(Cell_line == "GM24143") %>% group_by(PDALSeq_increase) %>% summarise(SUM = sum(sum))
df_sum_GM24143$Cell_line = "GM24143"
df_sum = bind_rows(df_sum_GM24385, df_sum_GM24143)
unique(df_sum$PDALSeq_increase)
df_sum$PDALSeq_increase = factor(df_sum$PDALSeq_increase, levels = c("S-phase", "G0/G1 arrest", "none"), labels = c("S-phase", "G0/G1 arrest", "neither"))


P_fraction = df_sum %>% ggplot(aes(y = Cell_line, x = SUM, fill = PDALSeq_increase)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = c("#E95562FF",  "#331068FF", "grey") )

PA = plot_grid(P, P, rel_widths = c(1.5,1), labels = c("A", "B"), label_size = 8)
P_out = plot_grid(PA, P_FC, P_fraction, ncol = 1, rel_heights = c(1,2.2,0.4), labels = c("","C", "D"), label_size = 8)

ggsave("revised_plots/Figure_S7_S-phase_contribution_to_PDAL-Seq_signal.svg", P_out, width = 6, height = 6.5, units = "in")
