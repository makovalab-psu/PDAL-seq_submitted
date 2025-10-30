library(car)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
#args = c("data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/1000000_nucleotides/Ptroglodytes_AG18354_PDALseq.bg", "data/16_lm_analysis_of_nonbmotifs/01_glm/Ptroglodytes_AG18354_PDALseq.csv")

Species = strsplit(strsplit(args[1], split = "/")[[1]][5], split = "_")[[1]][1]

####Read in data
df.bg = readr::read_delim(args[1], col_names = c("chr", "start", "end", "reads"))
df.bg$Norm_reads = df.bg$reads*10^6/sum(df.bg$reads)

NONB = c("APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z")

for (x in NONB){
	path_x = paste("data/04_Genome_wide_comparison_to_nonB_motifs/03_nonb_density/1000000_nucleotides/", Species, "/", x, ".bg", sep = "")
	df.x = readr::read_delim(path_x, col_names = c("chr", "start", "end", "count", "intersect", "length", "density"))
	df.bg$X = df.x$density
	colnames(df.bg)[which(colnames(df.bg) == "X")] = x
}

path_gc = paste("data/16_lm_analysis_of_nonbmotifs/00_genome_wide_GC_content/1000000_nucleotides/", Species, ".txt", sep = "")
df.GC = readr::read_delim(path_gc, col_names = c("GC_count"), delim = "\t")

df.bg$GC_content = df.GC$GC_count/(df.bg$end - df.bg$start + 1)

df = df.bg

# Null model
null_model <- glm(reads ~ 1, family = poisson(), data = df)

# nonB motif models
fit.APR = glm(reads ~ APR, family = poisson(), data = df)
fit.IR = glm(reads ~ IR, family = poisson(), data = df)
fit.CRU = glm(reads ~ CRU, family = poisson(), data = df)
fit.MR = glm(reads ~ MR, family = poisson(), data = df)
fit.TRI = glm(reads ~ TRI, family = poisson(), data = df)
fit.Z = glm(reads ~ Z, family = poisson(), data = df)
fit.STR = glm(reads ~ STR, family = poisson(), data = df) 
fit.DR = glm(reads ~ DR, family = poisson(), data = df)
fit.SLS = glm(reads ~ SLS, family = poisson(), data = df) 
fit.GQ = glm(reads ~ GQ, family = poisson(), data = df)
fit.GC_content = glm(reads ~ GC_content, family = poisson(), data = df)

# McFadden's R²
mcfadden_r2.APR <- 1 - (as.numeric(logLik(fit.APR)) / as.numeric(logLik(null_model)))
mcfadden_r2.IR <- 1 - (as.numeric(logLik(fit.IR)) / as.numeric(logLik(null_model)))
mcfadden_r2.CRU <- 1 - (as.numeric(logLik(fit.CRU)) / as.numeric(logLik(null_model)))
mcfadden_r2.MR <- 1 - (as.numeric(logLik(fit.MR)) / as.numeric(logLik(null_model)))
mcfadden_r2.TRI <- 1 - (as.numeric(logLik(fit.TRI)) / as.numeric(logLik(null_model)))
mcfadden_r2.Z <- 1 - (as.numeric(logLik(fit.Z)) / as.numeric(logLik(null_model)))
mcfadden_r2.STR <- 1 - (as.numeric(logLik(fit.STR)) / as.numeric(logLik(null_model)))
mcfadden_r2.DR <- 1 - (as.numeric(logLik(fit.DR)) / as.numeric(logLik(null_model)))
mcfadden_r2.SLS <- 1 - (as.numeric(logLik(fit.SLS)) / as.numeric(logLik(null_model)))
mcfadden_r2.GQ <- 1 - (as.numeric(logLik(fit.GQ)) / as.numeric(logLik(null_model)))
mcfadden_r2.GC_content <- 1 - (as.numeric(logLik(fit.GC_content)) / as.numeric(logLik(null_model)))

#Get and examine fit summaries
df.out = data.frame("Parameter" = c("APR","IR","CRU", "MR", "TRI", "Z", "STR", "DR", "SLS", "GQ", "GC_content"),
		     "mcfadden_r2" = c(mcfadden_r2.APR, mcfadden_r2.IR, mcfadden_r2.CRU, mcfadden_r2.MR, mcfadden_r2.TRI, mcfadden_r2.Z, mcfadden_r2.STR, mcfadden_r2.DR, mcfadden_r2.SLS, mcfadden_r2.GQ, mcfadden_r2.GC_content),
		     "Coef" =  c(coef(fit.APR)[2], coef(fit.IR)[2], coef(fit.CRU)[2], coef(fit.MR)[2], coef(fit.TRI)[2], coef(fit.Z)[2], coef(fit.STR)[2], coef(fit.DR)[2], coef(fit.SLS)[2], coef(fit.GQ)[2], coef(fit.GC_content)[2]))

write.csv(df.out, args[2], row.names = F)


