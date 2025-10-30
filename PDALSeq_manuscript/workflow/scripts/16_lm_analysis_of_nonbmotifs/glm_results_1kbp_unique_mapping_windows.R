library(dplyr)
library(viridis)
library(ggbeeswarm)
library(ggpubr)
library(cowplot)

#####Define variables####

Species_levels =  c("Hsapien","Ptroglodytes", "Ppaniscus",    "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus")
Species_labels = c("H. sapiens","P. troglodytes", "P. paniscus",    "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus")
Coef_levels = c("(Intercept)", "GC_content", "all", "APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS",  "GQ")
Coef_labels = c("Intercept", "GC_content", "all", "APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS",  "G4")
Fit_levels = c("all_GC", "all_no_GC", "dedup_GC", "dedup_no_GC")
Fit_labels = c("all with GC", "all no GC", "dedup with GC", "dedup no GC")

####Read in data####
df.ind = read.csv("results/16_lm_analysis_of_nonbmotifs/glm_unique_regions_1000_nucleotides_individual_parameters.csv")
df.global = read.csv("results/16_lm_analysis_of_nonbmotifs/glm_unique_regions_1000_nucleotides.csv")

####R2 data table and plot####

head(df.ind)
head(df.global)
df.R2 = df.global %>% filter(Fit == "all_GC") %>% filter(Coef == "(Intercept)") %>% select(Coef, mcfadden_r2, Experiment, Species, Bio_sample, Technique)
colnames(df.R2)[1] = "Parameter"
df.R2$Parameter = "all"
df.R2 = bind_rows(df.R2, df.ind %>% select(Parameter, mcfadden_r2, Experiment, Species, Bio_sample, Technique))

df.R2$Parameter = factor(df.R2$Parameter, levels = Coef_levels , labels = Coef_labels)
df.R2$Species = factor(df.R2$Species, levels = Species_levels, labels = Species_labels)

P.R2 = ggplot(data = df.R2, mapping = aes(x = Parameter, y = mcfadden_r2)) +
	facet_wrap(~Technique, nrow = 1) +
	geom_bar(stat = "summary", fun = mean, fill = "grey") +
	stat_summary(fun.data = mean_sd, geom = "errorbar", color = "black") +
	theme_classic() +
	scale_color_manual(values = viridis(length(unique(Species_levels)), option = "turbo")) +
	geom_hline(yintercept = 0) +
	theme(axis.title.x = element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1)) +
	      ylab("Mcfadden's R2")

df.R2.out = df.R2 %>% select(Species, Bio_sample, Technique, Parameter, mcfadden_r2) %>% arrange(Parameter) %>% arrange(Species) %>% arrange(Technique)

####VIF####

head(df.global)

df.global$Parameter = factor(df.global$Coef, levels = Coef_levels , labels = Coef_labels)
df.global$Species = factor(df.global$Species, levels = Species_levels, labels = Species_labels)
df.global$Fit = factor(df.global$Fit, levels = Fit_levels, labels = Fit_labels)

P.VIF = ggplot(data = df.global, mapping = aes(x = Parameter, y = vif, fill = Fit)) +
	facet_wrap(~Technique, nrow = 1) +
	geom_bar(stat = "summary", fun = mean, position = position_dodge()) +
	stat_summary(fun.data = mean_sd, geom = "errorbar", position = position_dodge()) +
	theme_classic() +
	scale_color_manual(values = viridis(length(unique(Species_levels)), option = "turbo")) +
	geom_hline(yintercept = 2) +
	theme(axis.title.x = element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1)) +
	      ylab("Varience inflation factor") +
	scale_fill_manual(values = viridis(4))
	
####Parameters####

P.Coef = ggplot(data = df.global %>% filter(Fit %in% c("dedup with GC", "dedup no GC")), mapping = aes(x = Parameter, y = Estimate, fill = Fit)) +
	facet_wrap(~Technique, nrow = 1) +
	geom_bar(stat = "summary", fun = mean, position = position_dodge()) +
	stat_summary(fun.data = mean_sd, geom = "errorbar", position = position_dodge()) +
	theme_classic() +
	scale_color_manual(values = viridis(length(unique(Species_levels)), option = "turbo")) +
	geom_hline(yintercept = 0) +
	theme(axis.title.x = element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1)) +
	      ylab("Fit coefficients (reads/kbp/motif)") +
	scale_fill_manual(values = viridis(4)[3:4])
	
###Compile plot####

P.out = plot_grid(P.R2, P.VIF, P.Coef, ncol = 1, labels = c("a", "b", "c"), align = "v")

ggsave("plots/Reads_vs_nonB_models/glm_poission_1kbp_windows_no_multimapping_regions/Summary_plot.svg", P.out, width = 6.5, height = 8, units = "in")

####Save tables####

write.csv(df.R2.out, "plots/Reads_vs_nonB_models/glm_poission_1kbp_windows_no_multimapping_regions/R2_by_coefficient.csv",  row.names = F)

df.global.r2 = df.global %>% filter(Coef == "(Intercept)") %>% select(Species, Bio_sample, Technique, Fit, mcfadden_r2) %>% arrange(Fit) %>% arrange(Species) %>% arrange(Technique)

write.csv(df.global.r2, "plots/Reads_vs_nonB_models/glm_poission_1kbp_windows_no_multimapping_regions/R2_by_global_fit.csv",  row.names = F)

write.csv(df.global %>% select(Species, Bio_sample, Technique, Fit, Parameter, Estimate,  Std..Error, z.value, Pr...z.., vif), "plots/Reads_vs_nonB_models/glm_poission_1kbp_windows_no_multimapping_regions/Global_fit_coefficients.csv",  row.names = F)

####Get R2 table for global fits####

df.global$Experiment_Fit = paste(df.global$Experiment, df.global$Fit)
l.df = df.global %>% group_by(Experiment_Fit) %>% group_split()
l.df.out = {}

for (i in 1:length(l.df)){
	l.df.out[[i]] = l.df[[i]][1,] %>% select(Species, Bio_sample, Fit, Technique, mcfadden_r2)
}

df.out = bind_rows(l.df.out) %>% arrange(Species) %>% arrange(Technique)

write.csv(df.out, "plots/Reads_vs_nonB_models/glm_poission_1kbp_windows_no_multimapping_regions/Global_fit_R2.csv",  row.names = F)

