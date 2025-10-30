library(car)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
#args = c("data/05_nonB_motif_cluster_analysis/11_compile_nonb_cluster_coverage/200/Hsapien_GM24385_TruSeq.csv")
df = readr::read_csv(args[1])
df$Length_N = df$end-df$start+1

# Null model
null_model <- glm(reads ~ 1, family = poisson(), data = df)

# nonB motif models
fit.Length_N = glm(reads ~ Length_N, family = poisson(), data = df)
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
mcfadden_r2.length_N <- 1 - (as.numeric(logLik(fit.Length_N)) / as.numeric(logLik(null_model)))
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
df.out = data.frame("Parameter" = c("Length_N", "APR","IR","CRU", "MR", "TRI", "Z", "STR", "DR", "SLS", "GQ", "GC_content"),
		     "mcfadden_r2" = c(mcfadden_r2.length_N, mcfadden_r2.APR, mcfadden_r2.IR, mcfadden_r2.CRU, mcfadden_r2.MR, mcfadden_r2.TRI, mcfadden_r2.Z, mcfadden_r2.STR, mcfadden_r2.DR, mcfadden_r2.SLS, mcfadden_r2.GQ, mcfadden_r2.GC_content),
		     "Coef" =  c(coef(fit.Length_N)[2], coef(fit.APR)[2], coef(fit.IR)[2], coef(fit.CRU)[2], coef(fit.MR)[2], coef(fit.TRI)[2], coef(fit.Z)[2], coef(fit.STR)[2], coef(fit.DR)[2], coef(fit.SLS)[2], coef(fit.GQ)[2], coef(fit.GC_content)[2]))


write.csv(df.out, args[2], row.names = F)
