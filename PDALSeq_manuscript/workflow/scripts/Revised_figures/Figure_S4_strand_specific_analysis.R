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

####Read and format data####

Species = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus")
Species_labels = c("human", "chimp", "bonobo", "gorilla", "B. orangutan", "S. orangutan", "siamang")
Species_palette = c("#9E9E9E", "#BE242F","#DC93B6", "#1C7D90", "#F2D983", "#FA9632", "#313757")
NONB_levels = rev(c("APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS",  "GQ", "all"))
NONB_labels = rev(c("APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS",  "G4", "all"))
NONB_keep = c("DR", "Z", "TRI", "CRU", "APR")

df = readr::read_csv("results/20_Strand_specific_analysis/Strand_specific_coverage.csv")
unique(df$Species)
df$Species = factor(df$Species, levels = Species, labels = Species_labels)

df$NonB_type = factor(df$NonB_type, levels = NONB_levels, labels = NONB_labels)
df_G4 = df  %>% filter(NonB_type %in% c("G4"))
df = df  %>% filter(NonB_type %in% NONB_keep)

df_G4 = df_G4 %>% select(Dataset, Species, Method, NonB_type, Position, Template_coverage, Nontemplate_coverage)
colnames(df_G4) = colnames(df)

df = bind_rows(df, df_G4)

unique(df$Dataset)
range(c(df$Forward_coverage,df$Reverse_coverage))

P_PDAL_F = ggplot(mapping = aes(x = Position, color = Species, group = Dataset)) +
  geom_vline(xintercept = 1) +
  geom_line(data = df %>% filter(Method == "PDALseq"), mapping = aes(y = Forward_coverage)) +
  facet_wrap(~NonB_type, nrow = 1) +
  scale_color_manual(values = Species_palette) +
  ylab("Mean coverage (CPM)") +
  xlab("Distance from start position (bp)")+
  ggtitle("PDAL-Seq forward") +
  scale_y_continuous(limits = c(0,0.14)) +
  scale_x_continuous(breaks = c(-500, 0, 500))

P_PDAL_R = ggplot(mapping = aes(x = Position, color = Species, group = Dataset)) +
  geom_vline(xintercept = 1) +
  geom_line(data = df %>% filter(Method == "PDALseq"), mapping = aes(y = Reverse_coverage)) +
  facet_wrap(~NonB_type, nrow = 1) +
  scale_color_manual(values = Species_palette) +
  ylab("Mean coverage (CPM)") +
  xlab("Distance from start position (bp)")+
  ggtitle("PDAL-Seq reverse") +
  scale_y_continuous(limits = c(0,0.14)) +
  scale_x_continuous(breaks = c(-500, 0, 500))


P_TRUE_F = ggplot(mapping = aes(x = Position, color = Species, group = Dataset)) +
  geom_vline(xintercept = 1) +
  geom_line(data = df %>% filter(Method != "PDALseq"), mapping = aes(y = Forward_coverage)) +
  facet_wrap(~NonB_type, nrow = 1) +
  scale_color_manual(values = Species_palette ) +
  ylab("Mean coverage (CPM)") +
  xlab("Distance from start position (bp)") +
  ggtitle("TruSeq forward") +
  scale_y_continuous(limits = c(0,0.14)) +
  scale_x_continuous(breaks = c(-500, 0, 500))

P_TRUE_R = ggplot(mapping = aes(x = Position, color = Species, group = Dataset)) +
  geom_vline(xintercept = 1) +
  geom_line(data = df %>% filter(Method != "PDALseq"), mapping = aes(y = Reverse_coverage)) +
  facet_wrap(~NonB_type, nrow = 1) +
  scale_color_manual(values = Species_palette) +
  ylab("Mean coverage (CPM)") +
  xlab("Distance from start position (bp)") +
  ggtitle("TruSeq reverse") +
  scale_y_continuous(limits = c(0,0.14)) +
  scale_x_continuous(breaks = c(-500, 0, 500))

P_PDAL = plot_grid(P_PDAL_F, P_PDAL_R, ncol = 1)
P_TRUE = plot_grid(P_TRUE_F, P_TRUE_R, ncol = 1)

P = plot_grid(P_PDAL, P_TRUE, ncol = 1, labels = c("A", "B"), label_size = 8)

ggsave("revised_plots/Figure_S4_Strand_specific_coverage.svg", P, width = 6.5, height = 5.5, units = "in")
