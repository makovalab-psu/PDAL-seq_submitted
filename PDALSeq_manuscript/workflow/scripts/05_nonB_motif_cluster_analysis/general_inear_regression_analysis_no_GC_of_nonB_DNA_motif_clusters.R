library(tidyverse)
library(car)

args <- commandArgs(trailingOnly = TRUE)

options(na.action = "na.fail")
args = "Hsapien_GM24385_PDALseq"
df = readr::read_csv("results/05_nonB_motif_cluster_analysis/nonb_cluster_coverage.csv") #%>% filter(Experiment == args[1])

P = ggplot(df, aes(x = read.count)) +
	geom_histogram() +
	theme_classic() +
	facet_wrap(~Experiment, scales = "free")
P
ggsave("plots/troubleshoot_glm/read_count_histogram.svg", P, width = 12.5, height = 6.5, units = "in")

P.log10 = ggplot(df, aes(x = log10(read.count))) +
	geom_histogram() +
	theme_classic() +
	facet_wrap(~Experiment, scales = "free")
P.log10
ggsave("plots/troubleshoot_glm/log10_read_count_histogram.svg", P.log10, width = 12.5, height = 6.5, units = "in")


for (i in 5:17){
	print(i)
	df.i = df[,c(i,22)]
	colnames(df.i) = c("x_data", "Experiment")
	x_lab = colnames(df)[i]
	log10_x_lab = paste("log10_", x_lab, sep = "")
	plot_path = paste("plots/troubleshoot_glm/", x_lab, ".svg", sep = "")
	log10_plot_path = paste("plots/troubleshoot_glm/", log10_x_lab, ".svg", sep = "")
	P = ggplot(df.i, aes(x = x_data)) +
	geom_histogram() +
	theme_classic() +
	facet_wrap(~Experiment, scales = "free")+
	xlab(x_lab)
	
	P.log10 = ggplot(df.i, aes(x = log10(x_data))) +
	geom_histogram() +
	theme_classic() +
	facet_wrap(~Experiment, scales = "free") +
	xlab(log10_x_lab)
	
	print(x_lab)
	ggsave(plot_path, P, width = 12.5, height = 6.5, units = "in")
	print(log10_x_lab)
	ggsave(log10_plot_path, P.log10, width = 12.5, height = 6.5, units = "in")
}

df %>% filter(Experiment == args[1])


#Normalize the data####

l.df = df %>%  group_by(Experiment) %>% group_split()

for (i in 1:length(l.df)){
	Norm_reads = (10^9)*l.df[[i]]$read.count/(sum(l.df[[i]]$read.count)*l.df[[i]]$length.N)
	l.df[[i]]$Norm_reads = Norm_reads
}

df = bind_rows(l.df)

P = ggplot(df, aes(x = Norm_reads)) +
	geom_histogram() +
	theme_classic() +
	facet_wrap(~Experiment, scales = "free")
ggsave("plots/troubleshoot_glm/Normalized_read_count_histogram.svg", P, width = 12.5, height = 6.5, units = "in")

P.nor_zoom = ggplot(df, aes(x = Norm_reads)) +
	geom_histogram() +
	theme_classic() +
	facet_wrap(~Experiment, scales = "free_y") +
	scale_x_continuous(limits = c(0,25))
	
ggsave("plots/troubleshoot_glm/Normalized_read_count_histogram.svg", P.nor_zoom, width = 12.5, height = 6.5, units = "in")

P.nor_zoom = ggplot(df, aes(x = Norm_reads)) +
	geom_histogram() +
	theme_classic() +
	facet_wrap(~Experiment, scales = "free_y") +
	scale_x_continuous(trans ="log10")
	
ggsave("plots/troubleshoot_glm/log10_Normalized_read_count_histogram.svg", P.nor_zoom, width = 12.5, height = 6.5, units = "in")

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

