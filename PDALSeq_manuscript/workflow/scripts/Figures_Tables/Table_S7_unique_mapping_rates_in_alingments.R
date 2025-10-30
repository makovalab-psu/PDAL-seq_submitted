library(tidyverse)

taxon = c("hominoidae", "hominidae", "homininae", "hominini")

l.df = {}
for (i in 1:length(taxon)){
	file_path = paste("data/07_PhyloHGMP_model/05_Unique_maps_inside_and_outside_alignment/", taxon[i], "/Unique_maps_inside_and_outside_alignment.csv", sep = "")
	df.i = read.csv(file_path)
	df.i$taxon = taxon[i]
	l.df[[i]] = df.i
}

df = bind_rows(l.df) %>% select(taxon, Species, Bio_sample, Technique, Unique.maps.alignment)

head(df)


SPECIES = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus")
SPECIES.label = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus")


BIO_SAMPLE = c("GM24385", "GM24143", "HEK293T", "CCL86-Raji", "AG18354", "AG18358", "KB8711", "KB14571", "KB3781", "F7549", "AG05252", "AG06213", "JAMBI", "Karenina")
BIO_SAMPLE.label = c("GM24385", "GM24143", "HEK293T", "CCL86-Raji", "AG18354", "AG18358", "KB8711", "KB14571", "KB3781", "F7549", "AG05252", "AG06213", "JAMBI", "Karenina")

df$Species = factor(df$Species, levels = SPECIES, labels = SPECIES.label) 
df$Bio_sample = factor(df$Bio_sample, levels = BIO_SAMPLE, labels = BIO_SAMPLE.label) 
df$taxon = factor(df$taxon, levels = taxon)

df.out = df %>% arrange(Bio_sample) %>% arrange(Species) %>% arrange(Technique) %>% arrange(taxon)

write.csv(df.out, file = "tables/Table_S7_Unique_maps_inside_alignment.csv", row.names = F)
