library(tidyverse)

df = read.csv("results/05_nonB_motif_cluster_analysis/linear_models.csv")

head(df)


SPECIES = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus")
SPECIES.label = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus")


BIO_SAMPLE = c("GM24385", "GM24143", "HEK293T", "CCL86-Raji", "AG18354", "AG18358", "KB8711", "KB14571", "KB3781", "F7549", "AG05252", "AG06213", "JAMBI", "Karenina")
BIO_SAMPLE.label = c("GM24385", "GM24143", "HEK293T", "CCL86-Raji", "AG18354", "AG18358", "KB8711", "KB14571", "KB3781", "F7549", "AG05252", "AG06213", "JAMBI", "Karenina")

df$Species = factor(df$Species, levels = SPECIES, labels = SPECIES.label) 
df$Bio_sample = factor(df$Bio_sample, levels = BIO_SAMPLE, labels = BIO_SAMPLE.label) 


df.out = df %>% arrange(Fit) %>% arrange(Bio_sample) %>% arrange(Species) %>% arrange(Technique) %>% select(Species, Bio_sample, Technique, Fit, X, Estimate, Std..Error, Pr...t.., VIF)

write.csv(df.out, file = "tables/Table_S4_Linear_model_parameters_and_test_for_multicolinearity.csv", row.names = F)
