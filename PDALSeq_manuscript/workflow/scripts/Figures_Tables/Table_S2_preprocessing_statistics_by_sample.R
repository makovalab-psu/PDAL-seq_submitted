library(tidyverse)

df = read.csv("results/02_preprocessing_statistics/Mapping_statistics_by_biological_sample.csv")

l.Experiment = strsplit(df$Experiment, split = "_")

v.species = c()
v.biosample = c()
v.technique = c()

for (i in 1:length(l.Experiment)){
	v.species[i] = l.Experiment[[i]][1]
	v.biosample[i] = l.Experiment[[i]][2]
	v.technique[i] = l.Experiment[[i]][3]
}

SPECIES = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus")
SPECIES.label = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus")
v.species = factor(v.species, levels = SPECIES, labels = SPECIES.label)

BIO_SAMPLE = c("GM24385", "GM24143", "HEK293T", "CCL86-Raji", "AG18354", "AG18358", "KB8711", "KB14571", "KB3781", "F7549", "AG05252", "AG06213", "JAMBI", "Karenina")
BIO_SAMPLE.label = c("GM24385", "GM24143", "HEK293T", "CCL86-Raji", "AG18354", "AG18358", "KB8711", "KB14571", "KB3781", "F7549", "AG05252", "AG06213", "JAMBI", "Karenina")
v.biosample = factor(v.biosample, levels = BIO_SAMPLE, labels = BIO_SAMPLE.label)

v.technique = factor(v.technique, levels = c("PDALseq", "TruSeq"), labels = c("PDAL-seq", "TruSeq"))

df$Species = v.species
df$Bio_sample = v.biosample
df$Technique = v.technique

df.out = df %>% arrange(Bio_sample) %>% arrange(Species) %>% arrange(Technique) %>% arrange(MapQ) %>% select(Species, Bio_sample, Technique, MapQ, mapped, secondary, uniquely.mapped)

df.out$Total = df.out$mapped + df.out$secondary
df.out$Fraction.unique = df.out$uniquely.mapped/df.out$Total
df.out$MapQ = factor(df.out$MapQ, levels = c("All_reads","High_MapQ"), labels = c("All reads", "High MapQ")) 


df.out %>% select(Species, Bio_sample, Technique, MapQ, mapped, secondary, Total, uniquely.mapped, Fraction.unique)




write.csv(df.out, file = "tables/Table_S2_preprocessing_statistics_by_sample.csv", row.names = F)
