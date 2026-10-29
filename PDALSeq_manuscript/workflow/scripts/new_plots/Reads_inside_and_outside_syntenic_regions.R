library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)
library(ggpubr)

theme_set(theme_classic(base_size = 6) +
            theme(text = element_text(color = "black", size = 6),
                  axis.text = element_text(color = "black", size = 6),
                  axis.title = element_text(color = "black", size = 6),
                  strip.text = element_text(color = "black", size = 6),
                  legend.text = element_text(color = "black", size = 6),
                  legend.title = element_text(color = "black", size = 6),
                  plot.title = element_text(color = "black", size = 6)))

list.files("results/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions")

df = readr::read_csv("results/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/PDAL_Seq_signal_in_and_out_of_syntenic_regions.csv")

head(df)
df$fraction_of_total_reads = df$Bases_in_region/(df$Total_bases)
df$Normalized_reads = df$Bases_in_region/(df$Total_bases*df$Region_nt)
unique(df$Region)
df$Region = factor(df$Region, levels = c("promoter", "utr5", "exon", "intron", "utr3", "censat", "nongenic_repeats", "nongenic_nonrepetitive"))

df$Method = NA
df$Method[which(grepl("PDALseq", df$Dataset))] = "PDAL-Seq"
df$Method[which(grepl("TruSeq", df$Dataset))] = "TruSeq"

Experiments = unique(df$Dataset)

l_df = {}
for (i in 1:length(Experiments)){
  l_df[[i]] = data.frame(df %>% filter(Dataset == Experiments[i]) %>% group_by(Syntenic) %>% summarise(sum = sum(fraction_of_total_reads)))
  l_df[[i]]$Experiment = Experiments[i]
  l_df[[i]]$Species = strsplit(Experiments[i], split = "_")[[1]][1] 
}
df_sum = bind_rows(l_df)
df_sum$Method = NA
df_sum$Method[which(grepl("PDALseq", df_sum$Experiment))] = "PDAL-Seq"
df_sum$Method[which(grepl("TruSeq", df_sum$Experiment))] = "TruSeq"
df_sum$Syntenic = factor(df_sum$Syntenic, levels = c("Syntenic", "Non-syntenic"))

df_sum$Species = factor(df_sum$Species, levels = c("Hsapien","Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus"), labels = c("H. sapiens","P. troglodytes", "P. paniscus", "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus"))
Species_palette = c("#9E9E9E", "#BE242F","#DC93B6", "#1C7D90", "#F2D983", "#FA9632", "#313757")

PA = ggplot(df_sum, aes(x = Syntenic, y = sum)) +
  geom_bar(stat = "summary", fun = "mean") +
  stat_summary(fun.data = mean_sd, geom = "errorbar", color = "black") +
  geom_jitter(mapping = aes(color = Species), size = 0.75) +
  facet_wrap(~Method, nrow = 1) + 
  scale_color_manual(values = Species_palette) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0,0.2,0.4,0.6, 0.8, 1.0)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  ylab("Fraction of reads")


df_GM24385 = df %>% filter(Dataset %in% c("Hsapien_GM24385_PDALseq", "Hsapien_GM24385_TruSeq"))
sum(df_GM24385 %>% filter(Dataset %in% c("Hsapien_GM24385_PDALseq")) %>% select(fraction_of_total_reads))
sum(df_GM24385 %>% filter(Dataset %in% c("Hsapien_GM24385_TruSeq")) %>% select(fraction_of_total_reads))

Region = unique(df$Region)
df_PDALSeq = df_GM24385 %>% filter(Method == "PDAL-Seq") %>% arrange(Syntenic) %>% arrange(Region)
df_TruSeq = df_GM24385 %>% filter(Method != "PDAL-Seq") %>% arrange(Syntenic) %>% arrange(Region)
Enrichment = df_PDALSeq$Normalized_reads/df_TruSeq$Normalized_reads
df_PDALSeq$Enrichment = Enrichment

df_PDALSeq$Syntenic = factor(df_PDALSeq$Syntenic, levels = c("Syntenic", "Non-syntenic"))

PB = ggplot(df_PDALSeq, aes(x = Region, y = log2(Enrichment), fill = Syntenic)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("black", "white")) +
  theme(axis.text.x = element_text(angle =90),
        axis.title.x = element_blank()) +
  ylab("log2(PDAL-Seq signal/TruSeq signal)") +
  ggtitle("GM24385-HG002")

P_fraction = ggplot(df, aes(x = Region, fill = Syntenic, y = fraction_of_total_reads)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge") +
  stat_summary(fun.data = mean_sd, geom = "errorbar", color = "black", position = "dodge") +
  facet_wrap(~Method, nrow = 1) + 
  scale_fill_manual(values = Species_palette) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  ylab("Fraction of reads")

P_norm = ggplot(df, aes(x = Region, fill = Syntenic, y = Normalized_reads)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge") +
  stat_summary(fun.data = mean_sd, geom = "errorbar", color = "black", position = "dodge") +
  facet_wrap(~Method, nrow = 1) + 
  scale_fill_manual(values = Species_palette) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  ylab("Fraction of reads")


ggsave("rebuttal_plots/Rebuttal_Figure_3_syntenic_nonsyntenic_analysis.svg", plot_grid(PA, PB, nrow = 1, labels = c("A", "B"), label_size = 6), width = 6, height = 3.5, units = "in")
