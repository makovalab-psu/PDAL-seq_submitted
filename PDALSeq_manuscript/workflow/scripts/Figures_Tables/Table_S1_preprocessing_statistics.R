library(tidyverse)

df = read.csv("results/02_preprocessing_statistics/Bam_flags_by_technical_replicate.csv")

SPECIES = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus")
SPECIES.label = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus")

df$Species = factor(df$Species, levels = SPECIES, labels = SPECIES.label)

BIO_SAMPLE = c("Hsapien_GM24385", "Hsapien_GM24143", "Hsapien_HEK293T", "Hsapien_CCL86-Raji", "Ptroglodytes_AG18354", "Ptroglodytes_AG18358", "Ppaniscus_KB8711", "Ppaniscus_KB14571", "Ggorilla_KB3781", "Ggorilla_F7549", "Ppygmaeus_AG05252", "Pabelii_AG06213", "Ssyndactylus_JAMBI", "Ssyndactylus_Karenina")
BIO_SAMPLE.label = c("GM24385", "GM24143", "HEK293T", "CCL86-Raji", "AG18354", "AG18358", "KB8711", "KB14571", "KB3781", "F7549", "AG05252", "AG06213", "JAMBI", "Karenina")

df$Bio_sample = factor(df$Bio_sample, levels = BIO_SAMPLE, labels = BIO_SAMPLE.label)

df.out = df %>% arrange(Bio_sample) %>% arrange(Technique) %>% select(Species, Bio_sample, Technique, Raw.fastq.million.reads, Trimmed.fastq.million.reads, Mapped.bam.primary.million.reads, Mapped.bam.supplementary.million.reads, Mapped.deduplicated.bam.primary.million.reads, Mapped.deduplicated.bam.supplementary.million.reads, Mapped.deduplicated.high.MapQ.bam.primary.million.reads, Mapped.deduplicated.high.MapQ.bam.supplementary.million.reads, Sequencing_sample) 

write.csv(df.out, file = "tables/Table_S1_preprocessing_statistics.csv", row.names = F)
