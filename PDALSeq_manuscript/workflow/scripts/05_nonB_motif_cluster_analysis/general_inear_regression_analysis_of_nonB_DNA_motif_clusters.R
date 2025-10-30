library(tidyverse)
library(car)

args <- commandArgs(trailingOnly = TRUE)

options(na.action = "na.fail")

df = readr::read_csv("results/05_nonB_motif_cluster_analysis/nonb_cluster_coverage.csv") %>% filter(Experiment == args[1])

#Normalize the data####

Norm_reads = (10^9)*df$read.count/sum(df$read.count)
df$Norm_reads = Norm_reads

#Fit the data to the global model
fit = glm(read.count ~ length.N + APR + IR + CRU + MR + TRI + Z + STR + DR + SLS + GQ + GC.content, df, family = poisson())

# Null model
null_model = glm(read.count ~ 1, family = poisson(), data = df)

# McFadden's R²
#R2=1−log-likelihood of the null model/log-likelihood of the fitted model​

mcfadden_r2 <- 1 - (as.numeric(logLik(fit)) / as.numeric(logLik(null_model)))
mcfadden_r2

#Test for multicolinearity with the variance inflation factor (VIF)
VIF = vif(fit)

df.out.1 = data.frame(coef(summary(fit)))
df.out.1$VIF = c(NA, VIF)
df.out.1$Experiment = args[1]
df.out.1$Fit = "Global"
df.out.1$X = rownames(df.out.1)

# Add Multiple R-squared and F-statistic as new columns
df.out.1$mcfadden_r2 = mcfadden_r2

#Fit the data to an optimized model
fit = glm(read.count ~ length.N + APR + CRU + TRI + Z + DR + GQ + GC.content, df, family = poisson())

#Test for multicolinearity with the variance inflation factor (VIF)
VIF = vif(fit)

# McFadden's R²
#R2=1−log-likelihood of the null model/log-likelihood of the fitted model​

mcfadden_r2 <- 1 - (as.numeric(logLik(fit)) / as.numeric(logLik(null_model)))
mcfadden_r2

df.out.2 = data.frame(coef(summary(fit)))
df.out.2$VIF = c(NA, VIF)
df.out.2$Experiment = args[1]
df.out.2$Fit = "Deduplicated parameters"
df.out.2$X = rownames(df.out.2)

# Add Multiple R-squared and F-statistic as new columns
df.out.2$mcfadden_r2= mcfadden_r2


#Consolidate

df.out = bind_rows(df.out.1, df.out.2)

####Write the averaged model coefficients####
write.csv(df.out, args[2], row.names = F)

