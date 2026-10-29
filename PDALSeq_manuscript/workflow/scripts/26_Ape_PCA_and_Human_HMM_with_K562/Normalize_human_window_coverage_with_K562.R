library(tidyverse)

existing_samples = c(
	"Hsapien_GM24385_PDALseq",
	"Hsapien_GM24143_PDALseq",
	"Hsapien_HEK293T_PDALseq",
	"Hsapien_CCL86-Raji_PDALseq"
)

normalize_bg = function(path){
	df = readr::read_delim(path, col_names = c("chr", "start", "end", "reads"),
	                       show_col_types = FALSE)
	Read_sum_billion = sum(df$reads, na.rm = TRUE) / 1e9
	Length_bp = df$end - df$start + 1
	rate = df$reads / Length_bp / Read_sum_billion
	Mean = mean(rate)
	SD   = mean(rate)
	(rate - Mean) / SD
}

l.Normalized = list()

for (i in seq_along(existing_samples)){
	path_i = paste0("data/09_Human_HMM/01_Remove_sex_chromosomes/", existing_samples[i], ".bg")
	l.Normalized[[i]] = normalize_bg(path_i)
}

l.Normalized[[5]] = normalize_bg(
	"data/26_Ape_PCA_and_Human_HMM_with_K562/04_remove_sex_chr/Hsapien_K562_PDALseq.bg"
)

df.out = data.frame(l.Normalized)
colnames(df.out) = c(existing_samples, "Hsapien_K562_PDALseq")

dir.create("data/26_Ape_PCA_and_Human_HMM_with_K562/05_normalized_coverage",
           recursive = TRUE, showWarnings = FALSE)
write.csv(df.out,
          "data/26_Ape_PCA_and_Human_HMM_with_K562/05_normalized_coverage/Normalized_human_window_coverage_with_K562.csv",
          row.names = FALSE)
