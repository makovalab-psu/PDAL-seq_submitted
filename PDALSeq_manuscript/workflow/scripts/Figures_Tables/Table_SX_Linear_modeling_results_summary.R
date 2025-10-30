library(tidyverse)

df = read.csv("results/05_nonB_motif_cluster_analysis/linear_models.csv")

head(df)


SPECIES = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus")
SPECIES.label = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus")


BIO_SAMPLE = c("GM24385", "GM24143", "HEK293T", "CCL86-Raji", "AG18354", "AG18358", "KB8711", "KB14571", "KB3781", "F7549", "AG05252", "AG06213", "JAMBI", "Karenina")
BIO_SAMPLE.label = c("GM24385", "GM24143", "HEK293T", "CCL86-Raji", "AG18354", "AG18358", "KB8711", "KB14571", "KB3781", "F7549", "AG05252", "AG06213", "JAMBI", "Karenina")

df$Species = factor(df$Species, levels = SPECIES, labels = SPECIES.label) 
df$Bio_sample = factor(df$Bio_sample, levels = BIO_SAMPLE, labels = BIO_SAMPLE.label) 


df.out = df %>% arrange(Fit) %>% arrange(Bio_sample) %>% arrange(Species) %>% arrange(Technique) %>% select(Species, Bio_sample, Technique, Fit, Multiple_R_squared, F_statistic, f_numdf, f_dendf, f_p_value)

df.out$Sample = paste(df.out$Species, df.out$Bio_sample, df.out$Technique, df.out$Fit, sep ="_")
df.summary = df.out %>% group_by(Sample) %>% summarize("Multiple_R_squared" = mean(Multiple_R_squared), "F_statistic" = mean(F_statistic), "f_numdf" = mean(f_numdf), "f_dendf" = mean(f_dendf), "f_p_value" = mean(f_p_value))

l.Sample = strsplit(df.summary$Sample, split = "_")
Species = c()
Bio_sample = c()
Technique = c()
Fit = c()

for (i in 1:nrow(df.summary)){
	Species[i] = l.Sample[[i]][1]
	Bio_sample[i] = l.Sample[[i]][2]
	Technique[i] = l.Sample[[i]][3]
	Fit[i] = l.Sample[[i]][4]
}

df.summary$Species = Species
df.summary$Bio_sample = Bio_sample
df.summary$Technique = Technique
df.summary$Fit = Fit
df.summary$Species = factor(df.summary$Species, levels = SPECIES.label, labels = SPECIES.label) 
df.summary$Bio_sample = factor(df.summary$Bio_sample, levels = BIO_SAMPLE, labels = BIO_SAMPLE.label) 
df.final = df.summary %>% arrange(Fit) %>% arrange(Bio_sample) %>% arrange(Species) %>% arrange(Technique) %>% select(Species, Bio_sample, Technique, Fit, Multiple_R_squared, F_statistic, f_numdf, f_dendf, f_p_value)

write.csv(df.final, file = "tables/Table_SX_Linear_model_statistics.csv", row.names = F)
