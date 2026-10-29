library(tidyverse)
library(viridis)

args <- commandArgs(trailingOnly = TRUE)

df = readr::read_csv(args[1], show_col_types = FALSE)

PDALseq = c(
	"Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq",
	"Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq",
	"Ptroglodytes_AG18354_PDALseq", "Ptroglodytes_AG18358_PDALseq",
	"Ppaniscus_KB8711_PDALseq", "Ppaniscus_KB14571_PDALseq",
	"Ggorilla_KB3781_PDALseq", "Ggorilla_F7549_PDALseq",
	"Ppygmaeus_AG05252_PDALseq", "Pabelii_AG06213_PDALseq",
	"Ssyndactylus_JAMBI_PDALseq", "Ssyndactylus_Karenina_PDALseq",
	"Hsapien_GM24385_TruSeq", "Ptroglodytes_AG18354_TruSeq",
	"Ppaniscus_KB8711_TruSeq", "Ggorilla_KB3781_TruSeq",
	"Ppygmaeus_AG05252_TruSeq", "Pabelii_AG06213_TruSeq",
	"Ssyndactylus_JAMBI_TruSeq",
	"Hsapien_K562_PDALseq"
)

df = df[, which(colnames(df) %in% PDALseq)]
df.clean = df[!rowSums(is.na(df)), ]

df.PCA = t(df.clean)

PCA = prcomp(df.PCA, scale = FALSE)

Experiments = colnames(df.clean)

var_explained = PCA$sdev^2 / sum(PCA$sdev^2)

P.scree = qplot(c(1:length(var_explained)), var_explained) +
	geom_line() +
	xlab("Principal Component") +
	ylab("Variance Explained") +
	ggtitle("Scree Plot")

l.Experiments = strsplit(Experiments, split = "_")

Species   = c()
Biosample = c()
Technique = c()

for (i in seq_along(l.Experiments)){
	Species[i]   = l.Experiments[[i]][1]
	Biosample[i] = l.Experiments[[i]][2]
	Technique[i] = l.Experiments[[i]][3]
}

df.out = data.frame(PCA$x)
df.out$Experiment = Experiments
df.out$Species    = factor(Species,
	levels = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus"),
	labels = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus"))
df.out$Biosample  = Biosample
df.out$Technique  = Technique

dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
write.csv(df.out, args[2], row.names = FALSE)
ggsave(args[3], P.scree)
