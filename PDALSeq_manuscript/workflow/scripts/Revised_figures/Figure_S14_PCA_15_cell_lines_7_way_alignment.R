library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)
library(ggrepel)

####Theme set####

theme_set(theme_classic(base_size = 6) +
            theme(text = element_text(color = "black", size = 6),
                  axis.text = element_text(color = "black", size = 6),
                  axis.title = element_text(color = "black", size = 6),
                  strip.text = element_text(color = "black", size = 6),
                  legend.text = element_text(color = "black", size = 6),
                  legend.title = element_text(color = "black", size = 6),
                  plot.title = element_text(color = "black", size = 6)))

df = readr::read_csv("results/26_Ape_PCA_and_Human_HMM_with_K562/hominoidae/PDALseq_PCA_in_alignment.csv")
colnames(df)
df$Species = factor(df$Species, levels = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus"), labels = c("human","chimpanzee", "bonobo", "gorilla", "B. orangutan", "S. orangutan", "siamang"))
Species_palette = c("#9E9E9E", "#BE242F","#DC93B6", "#1C7D90", "#F2D983", "#FA9632", "#313757")
unique(df$Biosample)
df$Cell_type = NA
df$Cell_type[which(df$Biosample %in% c("GM24385", "GM24143", "AG18354", "AG18358", "JAMBI", "Karenina"))] = "lymphoblast"
df$Cell_type[which(df$Biosample %in% c("KB8711", "KB14571", "KB3781", "F7549", "AG05252", "AG06213"))] = "fibroblast"
df$Cell_type[which(df$Biosample %in% c("HEK293T", "CCL86-Raji", "K562"))] = "other"
df$Cell_type = factor(df$Cell_type, levels = c("lymphoblast", "fibroblast", "other"))

P = df %>% ggplot(aes(x = PC1, y = PC2, color = Species, label = Biosample, shape = Cell_type)) +
  geom_point() +
  geom_text_repel(size = 2) +
  scale_color_manual(values = Species_palette)

ggsave("revised_plots/Figure_S10_PCA_15_cell_lines_7_way_alignment.svg", P, width = 4.5, height = 3, units = "in")
