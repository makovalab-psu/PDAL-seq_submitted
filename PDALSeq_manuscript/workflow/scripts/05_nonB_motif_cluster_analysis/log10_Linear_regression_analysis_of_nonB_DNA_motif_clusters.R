library(tidyverse)
library(car)

args <- commandArgs(trailingOnly = TRUE)

options(na.action = "na.fail")

df = readr::read_csv("results/05_nonB_motif_cluster_analysis/nonb_cluster_coverage.csv") %>% filter(Experiment == args[1])

#Normalize the data####

Norm_reads = (10^9)*df$read.count/(sum(df$read.count)*df$length.N)
df$log10_Norm_reads = log10(Norm_reads)
df = df[which(is.finite(df$log10_Norm_reads)),]

#Fit the data to the global model
fit = lm(log10_Norm_reads ~ APR + IR + CRU + MR + TRI + Z + STR + DR + SLS + GQ + GC.content, df)

#Test for multicolinearity with the variance inflation factor (VIF)
VIF = vif(fit)

df.out.1 = data.frame(coef(summary(fit)))
df.out.1$VIF = c(NA, VIF)
df.out.1$Experiment = args[1]
df.out.1$Fit = "Global"
df.out.1$X = rownames(df.out.1)

# Get the summary of the fitted model
fit_summary = summary(fit)

# Extract Multiple R-squared and F-statistic
multiple_r_squared = fit_summary$r.squared
f_statistic = fit_summary$fstatistic[1]  # F-statistic value
f_numdf = fit_summary$fstatistic[2]
f_dendf = fit_summary$fstatistic[3]

f_p_value = pf(f_statistic, f_numdf, f_dendf , lower.tail = FALSE)

# Add Multiple R-squared and F-statistic as new columns
df.out.1$Multiple_R_squared = multiple_r_squared
df.out.1$F_statistic = f_statistic
df.out.1$f_numdf = f_numdf
df.out.1$f_dendf = f_dendf 
df.out.1$f_p_value = f_p_value

#Fit the data to an optimized model
fit = lm(log10_Norm_reads ~ APR + CRU + TRI + Z + DR + GQ + GC.content, df)

#Test for multicolinearity with the variance inflation factor (VIF)
VIF = vif(fit)

df.out.2 = data.frame(coef(summary(fit)))
df.out.2$VIF = c(NA, VIF)
df.out.2$Experiment = args[1]
df.out.2$Fit = "Deduplicated parameters"
df.out.2$X = rownames(df.out.2)

# Get the summary of the fitted model
fit_summary = summary(fit)

# Extract Multiple R-squared and F-statistic
multiple_r_squared = fit_summary$r.squared
f_statistic = fit_summary$fstatistic[1]  # F-statistic value
f_numdf = fit_summary$fstatistic[2]
f_dendf = fit_summary$fstatistic[3]

f_p_value = pf(f_statistic, f_numdf, f_dendf , lower.tail = FALSE)

# Add Multiple R-squared and F-statistic as new columns
df.out.2$Multiple_R_squared = multiple_r_squared
df.out.2$F_statistic = f_statistic
df.out.2$f_numdf = f_numdf
df.out.2$f_dendf = f_dendf 

#Consolidate

df.out = bind_rows(df.out.1, df.out.2)

####Write the averaged model coefficients####
write.csv(df.out, args[2], row.names = F)

