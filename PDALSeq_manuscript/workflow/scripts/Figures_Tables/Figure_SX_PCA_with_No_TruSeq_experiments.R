library(dplyr)
library(ggplot2)
library(viridis)

args = commandArgs(trailingOnly = TRUE)

df = readr::read_csv("data/07_PhyloHGMP_model/alignment_coverage/hominoidae/alignment_coverage_matrix.csv")

head(df)
table(rowSums(is.na(df)))
nrow(df)

df = df[,-which(grepl("TruSeq", colnames(df)))]
df = df[,-which(colnames(df) %in% c("Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq"))]
df.PCA = t(df)

PCA = prcomp(df.PCA, scale = F)


#calculate total variance explained by each principal component
var_explained = PCA$sdev^2 / sum(PCA$sdev^2)

#create scree plot
P.scree = qplot(c(1:length(var_explained)), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") 

Experiments = colnames(df)

l.Experiments = strsplit(Experiments, split = "_")

Species = c()
Biosample = c()
Technique = c()

for (i in 1:length(l.Experiments)){
	Species[i] = l.Experiments[[i]][1]
	Biosample[i] = l.Experiments[[i]][2]
	Technique[i] = l.Experiments[[i]][3]
}

df.PCA = data.frame(PCA$x) 
df.PCA$Experiment = Experiments
df.PCA$Species = factor(Species, levels = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus"), labels = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus"))
df.PCA$Biosample = Biosample
df.PCA$Technique = Technique

P.PCA = ggplot(df.PCA, aes(x = PC1, y = PC2, color = Species, shape = Technique, label = Biosample)) +
	geom_point() +
	geom_text() +
	theme_classic() +
	scale_color_manual(values = viridis(length(unique(df.PCA$Species))+1, option = "turbo"))
	
ggsave("plots/Figure_SX_PCA_with_No_TruSeq_experiments/Figure_SX.svg", P.PCA, width = 5.5, height = 4, units = "in", scale = 1)

