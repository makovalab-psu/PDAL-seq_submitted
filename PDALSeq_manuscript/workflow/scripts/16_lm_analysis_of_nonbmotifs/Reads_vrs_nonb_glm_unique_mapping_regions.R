library(car)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
#args = c("data/16_lm_analysis_of_nonbmotifs/07_glm_matrix/1000_nucleotides/Hsapien_GM24385_PDALseq.csv", "GLM_test_analysis.csv")
df = readr::read_csv(args[1])

# Null model
null_model <- glm(reads ~ 1, family = poisson(), data = df)

# nonB motif models
fit.all.GC = glm(reads ~ APR + IR + CRU + MR + TRI + Z + STR + DR + SLS + GQ + GC_content, family = poisson(), data = df)
fit.all.no.GC = glm(reads ~ APR + IR + CRU + MR + TRI + Z + STR + DR + SLS + GQ, family = poisson(), data = df)
fit.dedup.GC = glm(reads ~ APR + CRU + TRI + Z + DR + GQ + GC_content, family = poisson(), data = df)
fit.dedup.no.GC = glm(reads ~ APR + CRU + TRI + Z + DR + GQ, family = poisson(), data = df)

#Get and examine fit summaries
Fit_summary.all.GC = summary(fit.all.GC)
Fit_summary.all.GC
Fit_summary.all.no.GC = summary(fit.all.no.GC)
Fit_summary.all.no.GC
Fit_summary.dedup.GC = summary(fit.dedup.GC)
Fit_summary.dedup.GC
Fit_summary.dedup.no.GC = summary(fit.dedup.no.GC)
Fit_summary.dedup.no.GC

# McFadden's R²
#R2=1−log-likelihood of the null model/log-likelihood of the fitted model​

mcfadden_r2.all.GC <- 1 - (as.numeric(logLik(fit.all.GC)) / as.numeric(logLik(null_model)))
mcfadden_r2.all.GC
mcfadden_r2.all.no.GC <- 1 - (as.numeric(logLik(fit.all.no.GC)) / as.numeric(logLik(null_model)))
mcfadden_r2.all.no.GC
mcfadden_r2.dedup.GC <- 1 - (as.numeric(logLik(fit.dedup.GC)) / as.numeric(logLik(null_model)))
mcfadden_r2.dedup.GC
mcfadden_r2.dedup.no.GC <- 1 - (as.numeric(logLik(fit.dedup.no.GC)) / as.numeric(logLik(null_model)))
mcfadden_r2.dedup.no.GC

#VIF 
vif.all.GC = vif(fit.all.GC)
vif.all.no.GC = vif(fit.all.no.GC)
vif.dedup.GC = vif(fit.dedup.GC)
vif.dedup.no.GC = vif(fit.dedup.no.GC)


#Write fit coefs, r2, vif

df.all.GC = data.frame(coef(Fit_summary.all.GC))
df.all.GC$mcfadden_r2 = mcfadden_r2.all.GC
df.all.GC$vif = c(NA, vif(fit.all.GC))
df.all.GC$Coef = row.names(df.all.GC)
df.all.GC$Fit = "all_GC"

df.all.no.GC = data.frame(coef(Fit_summary.all.no.GC))
df.all.no.GC$mcfadden_r2 = mcfadden_r2.all.no.GC
df.all.no.GC$vif = c(NA, vif(fit.all.no.GC))
df.all.no.GC$Coef = row.names(df.all.no.GC)
df.all.no.GC$Fit = "all_no_GC"

df.dedup.GC  = data.frame(coef(Fit_summary.dedup.GC))
df.dedup.GC$mcfadden_r2 = mcfadden_r2.dedup.GC
df.dedup.GC$vif = c(NA, vif(fit.dedup.GC))
df.dedup.GC$Coef = row.names(df.dedup.GC)
df.dedup.GC$Fit = "dedup_GC"

df.dedup.no.GC = data.frame(coef(Fit_summary.dedup.no.GC))
df.dedup.no.GC$mcfadden_r2 = mcfadden_r2.dedup.no.GC
df.dedup.no.GC$vif = c(NA, vif(fit.dedup.no.GC))
df.dedup.no.GC$Coef = row.names(df.dedup.no.GC)
df.dedup.no.GC$Fit = "dedup_no_GC"

df.out = bind_rows(df.all.GC, df.all.no.GC, df.dedup.GC, df.dedup.no.GC)

write.csv(df.out, args[2], row.names = F)


