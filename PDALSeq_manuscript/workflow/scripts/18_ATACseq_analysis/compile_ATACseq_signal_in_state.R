library(dplyr)

SRA_levels = c("SRR13051064", "SRR13051019", "SRR13051041",
               "SRR13051198", "SRR13051203", "SRR13051201")
Cell_line_levels = c("HG002-GM24385", "HG002-GM24385", "HG002-GM24385",
                     "HG004-GM24143", "HG004-GM24143", "HG004-GM24143")
Ape_states = 0:7

enrichment = function(Data_set = "SRR13051064", State = 3){
	File_path = paste("data/18_ATACseq_analysis/07_ape_state_ATACseq/state_", State, "/", Data_set, ".csv", sep = "")
	df = read.csv(File_path)
	Mean_observed = mean(df$Observed)
	Mean_null = mean(df$Null)
	Wilcox = wilcox.test(df$Observed, df$Null)
	Wilcox.p = Wilcox$p.value[[1]][1]
	df_out = data.frame(Data_set, State, Mean_observed, Mean_null, Wilcox.p)
	return(df_out)
}

enrichment_sra = function(sra = "SRR13051064"){
	l_df = {}
	for (i in 1:length(Ape_states)){
		l_df[[i]] = enrichment(sra, Ape_states[i])
	}
	df_out = bind_rows(l_df)
	return(df_out)
}

l_df = lapply(SRA_levels, enrichment_sra)

df = bind_rows(l_df)
df$Bio_sample = factor(df$Data_set, levels = SRA_levels, labels = Cell_line_levels)

df$Enrichment = df$Mean_observed/df$Mean_null
df$log2_Enrichment = log2(df$Enrichment)

write.csv(df, "results/18_ATACseq_analysis/ATACseq_read_enrichment.csv", row.names = F)

