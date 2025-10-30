
library(tidyverse)
library(viridis)

args <- commandArgs(trailingOnly = TRUE)

df = readr::read_csv(args[1])

if (args[4] == "hominoidae"){
	PDALseq = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq","Ptroglodytes_AG18354_PDALseq",  "Ptroglodytes_AG18358_PDALseq" , "Ppaniscus_KB8711_PDALseq", "Ppaniscus_KB14571_PDALseq", "Ggorilla_KB3781_PDALseq", "Ggorilla_F7549_PDALseq", "Ppygmaeus_AG05252_PDALseq", "Pabelii_AG06213_PDALseq", "Ssyndactylus_JAMBI_PDALseq", "Ssyndactylus_Karenina_PDALseq", "Hsapien_GM24385_TruSeq", "Ptroglodytes_AG18354_TruSeq", "Ppaniscus_KB8711_TruSeq", "Ggorilla_KB3781_TruSeq", "Ppygmaeus_AG05252_TruSeq", "Pabelii_AG06213_TruSeq", "Ssyndactylus_JAMBI_TruSeq")
}

if (args[4] == "hominidae"){
	PDALseq = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq","Ptroglodytes_AG18354_PDALseq",  "Ptroglodytes_AG18358_PDALseq" , "Ppaniscus_KB8711_PDALseq", "Ppaniscus_KB14571_PDALseq", "Ggorilla_KB3781_PDALseq", "Ggorilla_F7549_PDALseq", "Ppygmaeus_AG05252_PDALseq", "Pabelii_AG06213_PDALseq", "Hsapien_GM24385_TruSeq", "Ptroglodytes_AG18354_TruSeq", "Ppaniscus_KB8711_TruSeq", "Ggorilla_KB3781_TruSeq", "Ppygmaeus_AG05252_TruSeq", "Pabelii_AG06213_TruSeq")
}

if (args[4] == "homininae"){
	PDALseq = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq","Ptroglodytes_AG18354_PDALseq",  "Ptroglodytes_AG18358_PDALseq" , "Ppaniscus_KB8711_PDALseq", "Ppaniscus_KB14571_PDALseq", "Ggorilla_KB3781_PDALseq",  "Ggorilla_F7549_PDALseq", "Hsapien_GM24385_TruSeq", "Ptroglodytes_AG18354_TruSeq", "Ppaniscus_KB8711_TruSeq", "Ggorilla_KB3781_TruSeq")
}

if (args[4] == "hominini"){
	PDALseq = c("Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq","Ptroglodytes_AG18354_PDALseq",  "Ptroglodytes_AG18358_PDALseq" , "Ppaniscus_KB8711_PDALseq", "Ppaniscus_KB14571_PDALseq", "Hsapien_GM24385_TruSeq", "Ptroglodytes_AG18354_TruSeq", "Ppaniscus_KB8711_TruSeq")
}

df = df[,which(colnames(df) %in% PDALseq)]

nrow(df)
df.clean <- df[!rowSums(is.na(df)),]
nrow(df.clean)

df.PCA = t(df.clean)

PCA = prcomp(df.PCA, scale = F)

Experiments = colnames(df)
PCA$Experiment = Experiments

#calculate total variance explained by each principal component
var_explained = PCA$sdev^2 / sum(PCA$sdev^2)

#create scree plot
P.scree = qplot(c(1:length(var_explained)), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") 


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


write.csv(df.PCA , args[2], row.names = F)
ggsave(args[3], P.scree)

