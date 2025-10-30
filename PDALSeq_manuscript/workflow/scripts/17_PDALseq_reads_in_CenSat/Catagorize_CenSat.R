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

table(df_other_apes$Class)
df_other_apes$Class = factor(df_other_apes$Class, levels =  Ape_levels, labels = Ape_labels)
table(df_other_apes$Class)

df.all = bind_rows(df.human, df_other_apes)

df.all$Subclass[which(is.na(df.all$Subclass))] = "none"

options(scipen = 9999999)

for (x in c("Hsapien", Species)){
	for (y in levels(df.all$Class)){
		df.bed = df.all %>% filter(Species == x) %>% filter(Class == y)
		file_path = paste("data/17_PDALseq_reads_in_CenSat/01_censat_classes/", x, "/", y, ".bed", sep = "")
		write.table(df.bed, file_path, sep = "\t", col.names = F, row.names = F, quote = F)
	}
}

