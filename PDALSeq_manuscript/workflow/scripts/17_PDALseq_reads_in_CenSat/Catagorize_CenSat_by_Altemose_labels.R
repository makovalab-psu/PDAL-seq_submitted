library(dplyr)

read.chr.names = function(x = "resources/genomes/Ggorilla.chrom.sizes"){
	df = read.delim(x, header = F, sep = "\t")
	colnames(df) = c("chr", "end")
	return(df)
}

#The goal of this script is to make the satellite names in the censat files match Altemose science paper

Apes = c("resources/repeats/linneas_annotations/mGorGor1_v2.0_CenSat_v2.0.bed",
"resources/repeats/linneas_annotations/mPanPan1_v2.0_CenSat_v2.0.bed",
"resources/repeats/linneas_annotations/mPanTro3_v2.0_CenSat_v2.0.bed",
"resources/repeats/linneas_annotations/mPonAbe1_v2.0_CenSat_v2.0.bed",
"resources/repeats/linneas_annotations/mPonPyg2_v2.0_CenSat_v2.0.bed",
"resources/repeats/linneas_annotations/mSymSyn1_v2.1_CenSat_v2.0.bed")
Species = c("Ggorilla", "Ppaniscus", "Ptroglodytes", "Pabelii", "Ppygmaeus", "Ssyndactylus") 

####Human annotations from Altemose science####

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

human.active.HOR = c("S1C1/5/19H1L", "S2C2H1L", "S01/1C3H1L", "S2C4H1L", "S01C6H1L", "S1C7H1L", "S2C8H1L", "S2C9H1L", "S1C10H1L", "S3C11H1L", "S1C12H1L","S2C13/21H1L", "S2C14/22H1L", "S2C15H1L", "S1C16H1L", "S3C17H1L", "S2C18H1L", "S2C20H1L", "S3CXH1L", "S4CYH1L")

#Separate active HORs from other sats
df.human.active.HOR = df.human[which(df.human$Subclass %in% human.active.HOR),]
df.human.active.HOR$Altemose = "active_HOR"
df.human = df.human[-which(df.human$Subclass %in% human.active.HOR),]

#Define the character strings that we are using for all the other human classifications

human.inactive.HOR = c("hor")
human.d.m.HOR = c("dhor", "mon")
human.HSat1 = c("hsat1A", "hsat1B")
human.HSat2  = c("hsat2")
human.HSat3 = c("hsat3")
human.bSat = c("bsat")
human.Other.satellite = c("censat", "gsat")
human.rDNA = c("rDNA")
human.Transistion.region = c("ct")

#check that all results are covered

which(!df.human$Class %in% c(human.inactive.HOR, human.d.m.HOR, human.HSat1, human.HSat2, human.HSat3, human.bSat, human.Other.satellite, human.rDNA, human.rDNA, human.Transistion.region))

#Separate and label results

#"inactive_HOR"
df.human.inactive.HOR = df.human[which(df.human$Class %in% human.inactive.HOR),]
df.human.inactive.HOR$Altemose = "inactive_HOR"
#diverged_monomeric_HOR
df.human.d.m.HOR = df.human[which(df.human$Class %in% human.d.m.HOR),]
df.human.d.m.HOR$Altemose = "diverged_monomeric_HOR"
#HSat1
df.human.HSat1 = df.human[which(df.human$Class %in% human.HSat1),]
df.human.HSat1$Altemose = "HSat1"
#HSat2
df.human.HSat2 = df.human[which(df.human$Class %in% human.HSat2),]
df.human.HSat2$Altemose = "HSat2"
#HSat3
df.human.HSat3 = df.human[which(df.human$Class %in% human.HSat3),]
df.human.HSat3$Altemose = "HSat3"
#bSat = df.human[which(
df.human.bSat = df.human[which(df.human$Class %in% human.bSat),]
df.human.bSat$Altemose = "bSat"
#Other_satellite
df.human.Other.satellite = df.human[which(df.human$Class %in% human.Other.satellite),]
df.human.Other.satellite$Altemose = "Other_satellite"
#rDNA
df.human.rDNA = df.human[which(df.human$Class %in% human.rDNA),]
df.human.rDNA$Altemose = "rDNA"
#Transition_region
df.human.Transistion.region = df.human[which(df.human$Class %in% human.Transistion.region),]
df.human.Transistion.region$Altemose = "Transition_region"

df.human.out = bind_rows(df.human.active.HOR, df.human.inactive.HOR, df.human.d.m.HOR, df.human.HSat1, df.human.HSat2, df.human.HSat3, df.human.bSat, df.human.Other.satellite, df.human.rDNA, df.human.Transistion.region)

df.human.out$Species = "Hsapien"

#Make sure all the OG data is retained

df.test = read.delim("resources/repeats/linneas_annotations/chm13v2.0_censat_v2.0.bed", sep = "\t", skip = 1, header = F)[1:4]

nrow(df.human.out) == nrow(df.test)

write.csv(df.human.out, "results/17_PDALseq_reads_in_CenSat/Human_Cen_Sat_coordinates.csv", row.names = F)
####All others from Yoo nature


l.df_other_apes = {}

for (i in 1:length(Apes)){
	df = read.delim(Apes[i], header = F)[,1:4]
	colnames(df) = c("chr", "start", "end", "Satellite")
	chr_names_path = paste("resources/genomes/", Species[i], ".chrom.sizes", sep = "")
	df.chr = read.chr.names(chr_names_path)
	df = df[which(df$chr %in% df.chr$chr),]
	df$Species = Species[i]
	l.df_other_apes[[i]] = df
}

df_other_apes = bind_rows(l.df_other_apes)

Class = c()
Subclass = c()

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


#Define the character strings that we are using for all the other ape classifications

active.HOR = c("active_hor")
inactive.HOR = c("hor")
d.m.HOR = c("mon", "mon/hor", "dhor")
HSat1 = c("HSat1B", "HSat1A")
HSat2  = c("HSat2")
HSat3 = c("HSat3")
bSat = c("bSat")
Other.satellite = c("gSat", "cenSat")
rDNA = c("rDNA")
Transistion.region = c("ct")
#Stuff I'm not including
Misc = c("mixedSuperFamily", "subTerm", "GAP", "GAP,HSat3",        "HSat2_3", "subTerm_mon", "subTerm_mon/hor")

which(!df_other_apes$Class %in% c(active.HOR, inactive.HOR, d.m.HOR, HSat1, HSat2, HSat3, bSat, Other.satellite, rDNA, Transistion.region, Misc))

#"active_HOR"
df.active.HOR = df_other_apes[which(df_other_apes$Class %in% active.HOR),]
df.active.HOR$Altemose =  "active_HOR"
#"inactive_HOR"
df.inactive.HOR = df_other_apes[which(df_other_apes$Class %in% inactive.HOR),]
df.inactive.HOR$Altemose =  "inactive_HOR"
#diverged_monomeric_HOR
df.d.m.HOR = df_other_apes[which(df_other_apes$Class %in% d.m.HOR),]
df.d.m.HOR$Altemose =  "diverged_monomeric_HOR"
#HSat1
df.HSat1 = df_other_apes[which(df_other_apes$Class %in% HSat1),]
df.HSat1$Altemose =  "HSat1"
#HSat2
df.HSat2 = df_other_apes[which(df_other_apes$Class %in% HSat2),]
df.HSat2$Altemose =  "HSat2"
#HSat3
df.HSat3 = df_other_apes[which(df_other_apes$Class %in% HSat3),]
df.HSat3$Altemose =  "HSat3"
#bSat
df.bSat = df_other_apes[which(df_other_apes$Class %in% bSat),]
df.bSat$Altemose =  "bSat"
#Other_satellite
df.Other.satellite = df_other_apes[which(df_other_apes$Class %in% Other.satellite),]
df.Other.satellite$Altemose = "Other_satellite"
#rDNA
df.rDNA = df_other_apes[which(df_other_apes$Class %in% rDNA),]
df.rDNA$Altemose =  "rDNA"
#Transition_region
df.Transistion.region = df_other_apes[which(df_other_apes$Class %in% Transistion.region),]
df.Transistion.region$Altemose =  "Transition_region"
#Misc
df.Misc = df_other_apes[which(df_other_apes$Class %in% Misc),]
df.Misc$Altemose =  "Misc"

df.Apes = bind_rows(df.active.HOR, df.inactive.HOR, df.d.m.HOR, df.HSat1, df.HSat2, df.HSat3, df.bSat, df.Other.satellite, df.rDNA, df.Transistion.region, df.Misc)

#Check to make sure we did not loose data

nrow(df.Apes) == nrow(df_other_apes)

#Filter out Misc
df.Apes.out = df.Apes %>% filter(Altemose != "Misc")

#Filter filter out random utigs

df.Apes.out = df.Apes.out[which(!grepl("random", df.Apes.out$chr)),]

#Check to make sure all catagories in humans match all catagories in apes

which(!unique(df.Apes.out$Altemose) %in% unique(df.human.out$Altemose))
which(!unique(df.human.out$Altemose) %in% unique(df.Apes.out$Altemose))

df.all = bind_rows(df.human.out, df.Apes.out)

Altemose_class = unique(df.all$Altemose)

options(scipen = 9999999)

for (x in c("Hsapien", Species)){
	for (y in Altemose_class){
		df.bed = df.all %>% filter(Species == x) %>% filter(Altemose == y)
		file_path = paste("data/17_PDALseq_reads_in_CenSat/07_Altemose_paper_censat_classes/", x, "/", y, ".bed", sep = "")
		write.table(df.bed, file_path, sep = "\t", col.names = F, row.names = F, quote = F)
	}
}

