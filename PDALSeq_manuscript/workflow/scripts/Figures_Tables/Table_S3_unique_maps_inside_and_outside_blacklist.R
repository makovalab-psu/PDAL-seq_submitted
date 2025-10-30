library(tidyverse)

df = read.csv("results/03_Correlation_and_black_list_unmappable_regions/Unique_maps_inside_and_outside_blacklist.csv")

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
df$Fraction_genome_in_blacklist = df$Blacklist.length/df$Total.length

df.out = df %>% arrange(Bio_sample) %>% arrange(Species) %>% arrange(Technique) %>% select(Species, Bio_sample, Technique, Blacklist.length, Fraction_genome_in_blacklist, Unique.maps.blacklist, Unique.maps.outside.blacklist)


write.csv(df.out, file = "tables/Table_S3_Unique_maps_inside_and_outside_blacklist.csv", row.names = F)
