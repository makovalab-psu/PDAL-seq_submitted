
library(tidyverse)
library(viridis)

args <- commandArgs(trailingOnly = TRUE)

df = readr::read_csv(args[1])

if (args[3] == "hominoidae"){
	PDALseq = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq","Ptroglodytes_AG18354_PDALseq",  "Ptroglodytes_AG18358_PDALseq" , "Ppaniscus_KB8711_PDALseq", "Ppaniscus_KB14571_PDALseq", "Ggorilla_KB3781_PDALseq", "Ggorilla_F7549_PDALseq", "Ppygmaeus_AG05252_PDALseq", "Pabelii_AG06213_PDALseq", "Ssyndactylus_JAMBI_PDALseq", "Ssyndactylus_Karenina_PDALseq", "Hsapien_GM24385_TruSeq", "Ptroglodytes_AG18354_TruSeq", "Ppaniscus_KB8711_TruSeq", "Ggorilla_KB3781_TruSeq", "Ppygmaeus_AG05252_TruSeq", "Pabelii_AG06213_TruSeq", "Ssyndactylus_JAMBI_TruSeq")
}

if (args[3] == "hominidae"){
	PDALseq = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq","Ptroglodytes_AG18354_PDALseq",  "Ptroglodytes_AG18358_PDALseq" , "Ppaniscus_KB8711_PDALseq", "Ppaniscus_KB14571_PDALseq", "Ggorilla_KB3781_PDALseq", "Ggorilla_F7549_PDALseq", "Ppygmaeus_AG05252_PDALseq", "Pabelii_AG06213_PDALseq", "Hsapien_GM24385_TruSeq", "Ptroglodytes_AG18354_TruSeq", "Ppaniscus_KB8711_TruSeq", "Ggorilla_KB3781_TruSeq", "Ppygmaeus_AG05252_TruSeq", "Pabelii_AG06213_TruSeq")
}

if (args[3] == "homininae"){
	PDALseq = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq","Ptroglodytes_AG18354_PDALseq",  "Ptroglodytes_AG18358_PDALseq" , "Ppaniscus_KB8711_PDALseq", "Ppaniscus_KB14571_PDALseq", "Ggorilla_KB3781_PDALseq",  "Ggorilla_F7549_PDALseq", "Hsapien_GM24385_TruSeq", "Ptroglodytes_AG18354_TruSeq", "Ppaniscus_KB8711_TruSeq", "Ggorilla_KB3781_TruSeq")
}

if (args[3] == "hominini"){
	PDALseq = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq","Ptroglodytes_AG18354_PDALseq",  "Ptroglodytes_AG18358_PDALseq" , "Ppaniscus_KB8711_PDALseq", "Ppaniscus_KB14571_PDALseq", "Hsapien_GM24385_TruSeq", "Ptroglodytes_AG18354_TruSeq", "Ppaniscus_KB8711_TruSeq")
}

df = df[, which(colnames(df) %in% PDALseq)]

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
n_clusters = seq(2, 50, by = 4)

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

#Write the data
write.csv(df.scree, args[2], row.names = F)
