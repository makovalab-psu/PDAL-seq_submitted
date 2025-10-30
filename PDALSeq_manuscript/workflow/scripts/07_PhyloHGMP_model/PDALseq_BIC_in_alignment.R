
library(tidyverse)
library(viridis)
library(mclust)
library(cluster)

df = readr::read_csv("results/07_PhyloHGMP_model/alignment_coverage/hominoidea/alignment_coverage_matrix.csv")


nrow(df)
df.clean <- df[!rowSums(is.na(df)),]
nrow(df.clean)

####Test if possible####
#make this example reproducible
set.seed(1)

#perform k-means clustering with k = 4 clusters
km = kmeans(df.clean, centers = 4, nstart = 1)

####Scree plot####

# Decide how many clusters to look at
n_clusters = seq(2, 60, by = 5)

# Initialize total within sum of squares error: wss
wss = rep(NA, length(n_clusters))

set.seed(123)

# Look over 1 to n possible clusters
for (i in 1:length(n_clusters)) {
  print(i)
  # Fit the model: km.out
  km.out = kmeans(df.clean, centers = n_clusters[i], nstart = 20)
  # Save the within cluster sum of squares
  wss[i] = km.out$tot.withinss
}

# Produce a scree plot
df.scree = data.frame(n_clusters, wss)



#### 

d_clust = Mclust(df.clean, G=n_clusters)
