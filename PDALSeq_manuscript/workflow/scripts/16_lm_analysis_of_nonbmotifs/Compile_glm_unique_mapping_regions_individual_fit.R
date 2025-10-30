library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

Experiments = c('Ggorilla_F7549_PDALseq', 'Ppygmaeus_AG05252_TruSeq', 'Pabelii_AG06213_TruSeq', 'Hsapien_GM24143_PDALseq', 'Hsapien_GM24385_TruSeq', 'Hsapien_HEK293T_PDALseq', 'Pabelii_AG06213_PDALseq', 'Hsapien_CCL86-Raji_PDALseq', 'Ptroglodytes_AG18358_PDALseq', 'Ssyndactylus_JAMBI_TruSeq', 'Ppygmaeus_AG05252_PDALseq', 'Ppaniscus_KB14571_PDALseq', 'Ptroglodytes_AG18354_PDALseq', 'Ggorilla_KB3781_PDALseq', 'Ggorilla_KB3781_TruSeq', 'Ppaniscus_KB8711_PDALseq', 'Hsapien_GM24385_PDALseq', 'Ssyndactylus_JAMBI_PDALseq', 'Ptroglodytes_AG18354_TruSeq', 'Ssyndactylus_Karenina_PDALseq', 'Ppaniscus_KB8711_TruSeq')

l.df = {}

for (i in 1:length(Experiments)){
	v.experiment = strsplit(Experiments[i], split = "_")[[1]]
	Species = v.experiment[1]
	Bio_sample = v.experiment[2]
	Technique = v.experiment[3]
	file_path = paste("data/16_lm_analysis_of_nonbmotifs/08_glm_unique_regions_individual_parameters/1000_nucleotides/", Experiments[i], ".csv", sep = "")
	df.i = read.csv(file_path)
	df.i$Experiment = Experiments[i]
	df.i$Species = Species
	df.i$Bio_sample = Bio_sample
	df.i$Technique  = Technique
	l.df[[i]] = df.i
}

df.out = bind_rows(l.df)

write.csv(df.out, args[1], row.names = F)
