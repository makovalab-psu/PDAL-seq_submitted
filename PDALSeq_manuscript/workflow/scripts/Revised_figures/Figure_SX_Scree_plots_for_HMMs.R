
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

####Ape HMM scree plots####

list.files("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/")
df_all = readr::read_csv("data/07_PhyloHGMP_model/alignment_clustering/hominoidae/PDALseq_kmeans_scree_plot_in_alignment.csv")

P_all_Kmeans = df_all %>% ggplot(aes(x = n_clusters, y = wss)) +
  geom_point() +
  geom_line() +
  ggtitle("Ape HMM K-means scree plot") +
  xlab("Cluster number") +
  ylab("Within-Cluster Sum of Squares")


####Lympphoblast Kmeans####

df_lympho_Kmeans = readr::read_csv("data/29_Lymphoblast_specific_HMM/PDALseq_kmeans_scree_plot_in_alignment.csv")

P_lympho_Kmeans = df_lympho_Kmeans %>% ggplot(aes(x = n_clusters, y = wss)) +
  geom_point() +
  geom_line() +
  ggtitle("Lymphoblast HMM K-means scree plot") +
  xlab("Cluster number") +
  ylab("Within-Cluster Sum of Squares")

df_lympho_BIC = readr::read_csv("data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/BIC.csv")

P_lympho_BIC = df_lympho_BIC %>% ggplot(aes(x = States, y = BIC)) +
  geom_point() +
  geom_line() +
  ggtitle("Lymphoblast HMM BIC plot") +
  xlab("HMM state number") +
  ylab("Bayesian Information Criterion (BIC) ")

P_lympho = plot_grid(P_lympho_Kmeans, P_lympho_BIC, nrow = 1)

####Ape HMM - Fibroblast only####

df_fibro_Kmeans = readr::read_csv("data/28_Fibroblast_specific_HMM/PDALseq_kmeans_scree_plot_in_alignment.csv")

P_fibro_Kmeans = df_lympho_Kmeans %>% ggplot(aes(x = n_clusters, y = wss)) +
  geom_point() +
  geom_line() +
  ggtitle("Fibrooblast HMM K-means scree plot") +
  xlab("Cluster number") +
  ylab("Within-Cluster Sum of Squares")

df_fibro_BIC = readr::read_csv("data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/BIC.csv")

P_fibro_BIC = df_fibro_BIC %>% ggplot(aes(x = States, y = BIC)) +
  geom_point() +
  geom_line() +
  ggtitle("Fibroblast HMM BIC plot") +
  xlab("HMM state number") +
  ylab("Bayesian Information Criterion (BIC) ")

P_fibro = plot_grid(P_fibro_Kmeans, P_fibro_BIC, nrow = 1)

####Human HMM####

df_Human_Kmeans = readr::read_csv("results/26_Ape_PCA_and_Human_HMM_with_K562/Kmeans_scree_data.csv")

P_Human_Kmeans = df_Human_Kmeans %>% ggplot(aes(x = n_clusters, y = wss)) +
  geom_point() +
  geom_line() +
  ggtitle("Human HMM K-means scree plot") +
  xlab("Cluster number") +
  ylab("Within-Cluster Sum of Squares")

df_Human_BIC = readr::read_csv("data/26_Ape_PCA_and_Human_HMM_with_K562/Human_GHMM/BIC.csv")

P_Human_BIC = df_Human_BIC %>% ggplot(aes(x = States, y = BIC)) +
  geom_point() +
  geom_line() +
  ggtitle("Human HMM BIC plot") +
  xlab("HMM state number") +
  ylab("Bayesian Information Criterion (BIC) ")

P_Human = plot_grid(P_Human_Kmeans, P_Human_BIC, nrow = 1)



