library(dplyr)

Apes = c("resources/repeats/linneas_annotations/mGorGor1_v2.0_CenSat_v2.0.bed",
"resources/repeats/linneas_annotations/mPanPan1_v2.0_CenSat_v2.0.bed",
"resources/repeats/linneas_annotations/mPanTro3_v2.0_CenSat_v2.0.bed",
"resources/repeats/linneas_annotations/mPonAbe1_v2.0_CenSat_v2.0.bed",
"resources/repeats/linneas_annotations/mPonPyg2_v2.0_CenSat_v2.0.bed",
"resources/repeats/linneas_annotations/mSymSyn1_v2.1_CenSat_v2.0.bed")
Species = c("Ggorilla", "Ppaniscus", "Ptroglodytes", "Pabelii", "Ppygmaeus", "Ssyndactylus") 


df.human = read.delim("resources/repeats/linneas_annotations/chm13v2.0_censat_v2.0.bed", sep = "\t", skip = 1, header = F)[1:4]
colnames(df.human) = c("chr", "start", "end", "Satellite")
df.human$Species = "Hsapien"

Class = c()
Subclass = c()

for (i in 1:nrow(df.human)){
	v.Sat = strsplit(df.human$Satellite[i], split = '\\(')[[1]]
	if (length(v.Sat) > 1){
		Class[i] = strsplit(v.Sat[1], split = "_")[[1]][1]
		Subclass[i] = gsub("\\)", "", v.Sat[2])
	}else{
		Class[i] = strsplit(v.Sat[1], split = "_")[[1]][1]
		Subclass[i] = NA
	}
	print(paste(i, Class[i], Subclass[i], sep = "    "))
}

df.human$Class = Class
df.human$Subclass =  Subclass

l.df_other_apes = {}

for (i in 1:length(Apes)){
	df = read.delim(Apes[i], header = F)[,1:4]
	colnames(df) = c("chr", "start", "end", "Satellite")
	df$Species = Species[i]
	l.df_other_apes[[i]] = df
}

df_other_apes = bind_rows(l.df_other_apes)

for (i in 1:nrow(df_other_apes)){
	v.Sat = strsplit(df_other_apes$Satellite[i], split = '\\(')[[1]]
	if (length(v.Sat) > 1){
		Class[i] = v.Sat[1]
		Subclass[i] = gsub("\\)", "", v.Sat[2])
	}else{
		Class[i] = v.Sat[1]
		Subclass[i] = NA
	}
	print(paste(i, Class[i], Subclass[i], sep = "    "))
}

df_other_apes$Class = Class
df_other_apes$Subclass =  Subclass

unique(df.human$Class)[order(unique(df.human$Class))]
unique(df_other_apes$Class)[order(unique(df_other_apes$Class))]

Human_levels  = c("active_hor", "bsat", "censat", "ct", "dhor", "GAP", "GAP,HSat3", "gsat", "hor", "hsat1A", "hsat1B",  "hsat2", "HSat2_3", "hsat3", "mixedSuperFamily", "mon", "mon/hor", "rDNA", "subTerm", "subTerm_mon", "subTerm_mon/hor")
Ape_levels = c("active_hor", "bSat", "cenSat", "ct", "dhor", "GAP", "GAP,HSat3", "gSat", "hor", "HSat1A", "HSat1B", "HSat2", "HSat2_3", "HSat3", "mixedSuperFamily", "mon", "mon/hor", "rDNA", "subTerm", "subTerm_mon", "subTerm_mon/hor")
Ape_labels = c("active_hor", "bSat", "cenSat", "ct", "dhor", "GAP", "GAP,HSat3", "gSat", "hor", "HSat1A", "HSat1B", "HSat2", "HSat2_3", "HSat3", "mixedSuperFamily", "mon", "mon-hor", "rDNA", "subTerm", "subTerm_mon", "subTerm_mon-hor")


table(df.human$Class)
df.human$Class = factor(df.human$Class, levels = Human_levels, labels =  Ape_labels)
table(df.human$Class)

df.all = bind_rows(df.human, df_other_apes)

df.all$Subclass[which(is.na(df.all$Subclass))] = "none"

####Count subclasses in classes####

l.df.all = df.all %>% group_by(Class) %>% group_split()
l.df.count = {}
for (i in 1:length(l.df.all)){
	df.i = l.df.all[[i]] %>% group_by(Subclass)  %>% summarize(BP = sum(end - start + 1))
	df.i$Class = l.df.all[[i]]$Class[1]
	l.df.count[[i]] = df.i
}

df.count = bind_rows(l.df.count)

library(ggplot2)

ggplot(df.count, aes(x = BP)) +
	geom_histogram() +
	theme_classic() +
	facet_wrap(~Class, scales = "free")
	
data.frame(table(df.count$Class)) %>% arrange(Freq)
Excluded_class = c("GAP", "GAP,HSat3", "HSat1A", "HSat2_3", "mon", "rDNA", "subTerm_mon", "subTerm_mon/hor", "ct")
df.subclass = df.count %>% filter(!Class %in% Excluded_class) %>% filter(BP >= 10^6) %>% arrange(-BP)

#Clean up so the files write nice

df.subclass$Subclass = gsub("/", "-", df.subclass$Subclass)
df.subclass$Class = gsub("/", "-", df.subclass$Class)

df.subclass$file_path = paste(df.subclass$Class, df.subclass$Subclass, sep = "/")

write.csv(df.subclass, "data/17_PDALseq_reads_in_CenSat/00_scratch/Filtered_subclasses_for_analysis.csv", row.names = F)

