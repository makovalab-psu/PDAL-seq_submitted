library(ggplot2)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)

df = read.csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/Transition_probabilities.csv")
df.arrange = readr::read_csv("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/Average_PDALseq_in_GHMM_states.csv") %>% arrange(MEAN)
States = paste("State_", unique(df.arrange$state), sep = "")


df$from_state = factor(as.character(df$from_state), levels = States, labels = gsub("State_", "", States))
df$to_state = factor(as.character(df$to_state), levels = States, labels = gsub("State_", "", States))
P.out = ggplot(df, aes(x = from_state, y = to_state, fill = transition_probability)) +
	geom_tile() +
	scale_fill_gradient(low = "white", high = "red", name = "transition probability") +
	theme_classic() +
	xlab("from state") +
	ylab("to state")

ggsave(args[1], P.out, height = 5, width = 6, scale = 1, limitsize = FALSE)
