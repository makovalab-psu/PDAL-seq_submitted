library(dplyr)
library(ggplot2)
library(cowplot)

df = readr::read_csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/GHMM_states.csv")

df$hidden_state = factor(as.character(df$hidden_state), levels = c("2", "5", "0", "6", "1", "3", "7", "4"))

Rows = 10000:10300

P.human = ggplot(df[Rows,], aes(x = Rows, fill = hidden_state)) +
	geom_bar(mapping = aes(y = Hsapien_GM24385_PDALseq), stat = "identity") +
	scale_fill_manual(values = c("#386350", "#386350", "#386350", "#B6DBFF", "#B6DBFF", "#882D71", "#DB5829", "#F4A637")) +
	theme_classic()+
	theme(axis.title.y = element_blank(),
	      axis.title.x = element_blank()) +
	ggtitle("human")
	
P.chimp = ggplot(df[Rows,], aes(x = Rows, fill = hidden_state)) +
	geom_bar(mapping = aes(y = Ptroglodytes_AG18354_PDALseq), stat = "identity") +
	scale_fill_manual(values = c("#386350", "#386350", "#386350", "#B6DBFF", "#B6DBFF", "#882D71", "#DB5829", "#F4A637")) +
	theme_classic()+
	theme(axis.title.y = element_blank(),
	      axis.title.x = element_blank()) +
	ggtitle("chimpanzee")
		
P.Borang = ggplot(df[Rows,], aes(x = Rows, fill = hidden_state)) +
	geom_bar(mapping = aes(y = Ppygmaeus_AG05252_PDALseq), stat = "identity") +
	scale_fill_manual(values = c("#386350", "#386350", "#386350", "#B6DBFF", "#B6DBFF", "#882D71", "#DB5829", "#F4A637")) +
	theme_classic() +
	theme(axis.title.y = element_blank(),
	      axis.title.x = element_blank()) +
	ggtitle("orangutan")
	
P = plot_grid(P.human, P.chimp, P.Borang, ncol = 1, align = "h")

ggsave("new_plots/HMM_plot.svg", P, width = 5, height = 5, units = "in")
