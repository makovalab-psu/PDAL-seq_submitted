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

add_phase_cell_line = function(df = df_nonB_cov){
  Phase = c()
  Cell_line = c()
  l_ann = strsplit(df$ann, split = "_")
  for (i in 1:length(l_ann)){
    Phase[i] = l_ann[[i]][1]
    Cell_line[i] = l_ann[[i]][2]
  }
  df$Phase = factor(Phase, levels = c("Sphase", "G0G1", "none"), labels = c("S-phase", "G0/G1 arrest", "neither"))
  df$Cell_line = factor(Cell_line, levels = c("GM24385", "GM24143"), labels = c("GM24385-HG002", "GM24143-HG004"))
  return(df)
}

NonB_levels = c("GQ", "SLS", "DR", "STR", "Z", "TRI", "MR", "CRU", "IR", "APR")
NonB_labels = c("G4", "SLS", "DR", "STR", "Z", "TRI", "MR", "CRU", "IR", "APR")

Fun_annotation_levels = c("promoter", "utr5", "exon", "intron", "utr3", "coreorigins", "stochorigins", "censat", "repeats")
Fun_annotation_labels = c("promoter", "5'-UTR", "exon", "intron", "3'-UTR", "core replication origins", "stochastic replication origins", "censat", "repeats")

Phase_pallette = c("#E95562FF",  "#331068FF", "grey")

####Read in data####

df_cov = add_phase_cell_line(readr::read_csv("results/25_Sphase_contribution/Sphase_nonb_func_coverage.csv"))
df_enr = add_phase_cell_line(readr::read_csv("results/25_Sphase_contribution/Sphase_nonb_func_enrichment.csv"))

####Non-B DNA element enrichment####

df_nonB_cov = df_cov %>% filter(feat %in% NonB_levels)
df_nonB_cov$feat = factor(df_nonB_cov$feat, levels = rev(NonB_levels), labels = rev(NonB_labels))

P_nonB_fraction = df_nonB_cov %>% ggplot(aes(x = 100*fraction_covered, y = feat, fill = Phase)) +
  geom_bar(stat = "summary", fun = mean, position = "dodge") +
  stat_summary(fun.data = mean_sd, geom = "errorbar", color = "black", position = "dodge") +
  scale_fill_manual(values = Phase_pallette, name = "PDAL-Seq increases in:") +
  ylab("non-B DNA motif") +
  xlab("Fraction of cell cycle specific region intersecting with non-B motifs (%)")  +
  theme(legend.position = "bottom") +
  ggtitle("non-B motif coverage") +
  scale_x_continuous(limits = c(0, 75))

df_nonB_enr = df_enr %>% filter(feat %in% NonB_levels)
df_nonB_enr$feat = factor(df_nonB_enr$feat, levels = rev(NonB_levels), labels = rev(NonB_labels))

P_nonB_enr = df_nonB_enr %>% ggplot(aes(x = Enrichment, y = feat, fill = Phase)) +
  geom_bar(stat = "summary", fun = mean, position = "dodge") +
  stat_summary(fun.data = mean_sd, geom = "errorbar", color = "black", position = "dodge") +
  scale_fill_manual(values = Phase_pallette, name = "PDAL-Seq increases in:") +
  ylab("non-B DNA motif") +
  xlab("Enrichment of non-B motifs in cell cycle specific region") +
  geom_vline(xintercept = 1) +
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank())  +
  ggtitle("non-B motif enrichment") +
  scale_x_continuous(limits = c(0, 5))

####Functional element enrichment####

df_fun_cov = df_cov %>% filter(feat %in% Fun_annotation_levels)
df_fun_cov$feat = factor(df_fun_cov$feat, levels = rev(Fun_annotation_levels), labels = rev(Fun_annotation_labels))

P_fun_fraction = df_fun_cov %>% ggplot(aes(x = 100*fraction_covered, y = feat, fill = Phase)) +
  geom_bar(stat = "summary", fun = mean, position = "dodge") +
  stat_summary(fun.data = mean_sd, geom = "errorbar", color = "black", position = "dodge") +
  scale_fill_manual(values = Phase_pallette, name = "PDAL-Seq increases in:") +
  ylab("Functional element") +
  xlab("Fraction of cell cycle specific region intersecting with functional elements (%)")  +
  theme(legend.position = "bottom") +
  ggtitle("Functional element coverage") +
  scale_x_continuous(limits = c(0, 75))

df_fun_enr = df_enr %>% filter(feat %in% Fun_annotation_levels)
df_fun_enr$feat = factor(df_fun_enr$feat, levels = rev(Fun_annotation_levels), labels = rev(Fun_annotation_labels))

P_fun_enr = df_fun_enr %>% ggplot(aes(x = Enrichment, y = feat, fill = Phase)) +
  geom_bar(stat = "summary", fun = mean, position = "dodge") +
  stat_summary(fun.data = mean_sd, geom = "errorbar", color = "black", position = "dodge") +
  scale_fill_manual(values = Phase_pallette, name = "PDAL-Seq increases in:") +
  ylab("Functional element") +
  xlab("Enrichment of Functional elements in cell cycle specific region") +
  geom_vline(xintercept = 1) +
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank())  +
  ggtitle("Functional element enrichment") +
  scale_x_continuous(limits = c(0, 5))

####Ape HMM state enrichment####

df_HMM_cov = add_phase_cell_line(readr::read_csv("results/25_Sphase_contribution/Cell_cycle_sensitive_region_coverage_of_Ape_HMM_states.csv"))
df_HMM_enr = add_phase_cell_line(readr::read_csv("results/25_Sphase_contribution/Cell_cycle_sensitive_region_enrichment_in_Ape_HMM_states.csv"))
State.levels = c("2", "5", "0", "6", "1", "3", "7", "4") 
State.labels = rev(c("8-origin", "7-recomb", "6-gene", "5-variable", "4-variable", "3-low", "2-low", "1-low"))

df_HMM_cov$State = factor(df_HMM_cov$State, levels = State.levels, labels = State.labels)
df_HMM_enr$State = factor(df_HMM_enr$State, levels = State.levels, labels = State.labels)

P_HMM_fraction = df_HMM_cov %>% ggplot(aes(x = 100*fraction_covered, y = State, fill = Phase)) +
  geom_bar(stat = "summary", fun = mean, position = "dodge") +
  stat_summary(fun.data = mean_sd, geom = "errorbar", color = "black", position = "dodge") +
  scale_fill_manual(values = Phase_pallette, name = "PDAL-Seq increases in:") +
  ylab("Ape HMM state") +
  xlab("Fraction of HMM state intersecting with cell cycle specific region (%)")  +
  theme(legend.position = "bottom") +
  ggtitle("Cell cycle specific region coverage") 

df_fun_enr = df_enr %>% filter(feat %in% Fun_annotation_levels)
df_fun_enr$feat = factor(df_fun_enr$feat, levels = rev(Fun_annotation_levels), labels = rev(Fun_annotation_labels))

P_HMM_enr = df_HMM_enr %>% ggplot(aes(x = Enrichment, y = State, fill = Phase)) +
  geom_bar(stat = "summary", fun = mean, position = "dodge") +
  stat_summary(fun.data = mean_sd, geom = "errorbar", color = "black", position = "dodge") +
  scale_fill_manual(values = Phase_pallette, name = "PDAL-Seq increases in:") +
  ylab("Functional element") +
  xlab("Enrichment of cell cycle specific region in HMM states") +
  geom_vline(xintercept = 1) +
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank())  +
  ggtitle("Cell cycle specific region enrichment") 

####Assemble plots####

PA = plot_grid(P_nonB_fraction, P_nonB_enr, P_fun_fraction, P_fun_enr, align = "v")
PA

PB = plot_grid(P_HMM_fraction, P_HMM_enr, align = "h")

P = plot_grid(PA, PB, ncol = 1, labels = c("A", "B"), label_size = 10, rel_heights = c(2:1))

ggsave("revised_plots/Figure_S8_functional_associations_with_cell_cycle_specific_PDAL-Seq_regions.svg", P, width = 6.5, height = 6.5, units = "in")
