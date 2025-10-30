library(tidyverse)
library(viridis)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)

#args = c("data/07_PhyloHGMP_model/alignment_coverage/homininae/alignment_coverage_matrix.csv", "data/07_PhyloHGMP_model/alignment_clustering/homininae/PDALseq_correllation_in_alignment.csv", "homininae")

df = readr::read_csv(args[1])

head(df)

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


Experiments = colnames(df)

df.cor = data.frame(expand.grid(PDALseq, PDALseq))
colnames(df.cor) = c("Sample_1", "Sample_2")
pearson = c()
spearman = c()
for (i in 1:nrow(df.cor)){
	print(i)
	x = which(Experiments == df.cor$Sample_1[i])
	y = which(Experiments == df.cor$Sample_2[i])
	t.pearson = cor.test(df[[x]], df[[y]], method = "pearson")
	t.spearman = cor.test(df[[x]], df[[y]], method = "spearman")
	pearson[i] = t.pearson$estimate 
	spearman[i] = t.spearman$estimate
}

df.cor$pearson = pearson
df.cor$spearman = spearman

write.csv(df.cor, args[2], row.names = F)

df.cor$Sample_1 = factor(df.cor$Sample_1, levels = PDALseq)
df.cor$Sample_2 = factor(df.cor$Sample_2, levels = PDALseq)

P.pearson = ggplot(df.cor, aes(x = Sample_1, y = Sample_2, fill = pearson)) +
	geom_tile() +
	scale_fill_viridis(option = "turbo", limits = c(-1, 1)) +
	theme_classic() +
	coord_fixed() +
	theme(axis.text.y = element_blank(),
	      axis.text.x = element_text(angle = 90))
	      
P.spearman = ggplot(df.cor, aes(x = Sample_1, y = Sample_2, fill = spearman)) +
	geom_tile() +
	scale_fill_viridis(option = "turbo", limits = c(-1, 1)) +
	theme_classic() +
	coord_fixed() +
	theme(axis.text.y = element_blank(),
	      axis.text.x = element_text(angle = 90))
  
P = plot_grid(P.pearson, P.spearman, nrow = 1)

ggsave(gsub(".csv", ".png", args[2]), P, width = 5, height = 4, scale = 2)
