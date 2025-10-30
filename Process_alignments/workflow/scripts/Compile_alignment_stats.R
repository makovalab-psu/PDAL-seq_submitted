library(dplyr)

Taxon = c("hominoidae", "hominidae", "homininae", "hominini")
Species = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus")

df = data.frame(expand.grid(Taxon, Species))
colnames(df) = c("Taxon", "Species")
head(df)

Blocks = c()
Merged_blocks = c()
BP = c()
Merged_BP = c()

for (i in 1:nrow(df)){
  file_path = paste("data/check_alignment_stats/", df$Taxon[i], "/", df$Species[i], ".txt", sep = "")
  CON = file(file_path)
  Lines = readLines(CON)
  close(CON)
  l.Lines = strsplit(Lines, split = " ")
  Blocks[i] = as.numeric(l.Lines[[1]][3])
  Merged_blocks[i] = as.numeric(l.Lines[[3]][3])
  BP[i] = as.numeric(l.Lines[[2]][2])
  Merged_BP[i] = as.numeric(l.Lines[[4]][2])
}

df$Blocks = Blocks
df$Merged_blocks = Merged_blocks
df$BP = BP
df$Merged_BP = Merged_BP

head(df)

write.csv(df, "data/check_alignment_stats/summary.csv", row.names = F)

