

#Low complexity
#"Low_complexity"
cat data/12_PDALseq/rm_annotation/Low_complexity.bed | bedtools sort > data/14_Repeat_masker_class_enrichment/01_compile_annotations/Low_complexity.bed

#Simple repeat
#"Simple_repeat"
cat data/12_PDALseq/rm_annotation/Simple_repeat.bed | bedtools sort > data/14_Repeat_masker_class_enrichment/01_compile_annotations/Simple_repeat.bed

#Satelitte
#"Satellite"
#"Satellite/acro"
#"Satellite/centr"
#"Satellite/subtelo"
cat data/12_PDALseq/rm_annotation/Satellite/*.bed data/12_PDALseq/rm_annotation/Satellite.bed | bedtools sort > data/14_Repeat_masker_class_enrichment/01_compile_annotations/Satelitte.bed

#DNA
#1"DNA"
#2"DNA?"
#3"DNA/Crypton"
#4"DNA/Crypton-A"
#5"DNA/hAT"
#6"DNA/hAT?"
#7"DNA/hAT-Ac"
#8"DNA/hAT-Blackjack"
#9"DNA/hAT-Charlie"
#10"DNA/hAT-hAT19"
#11"DNA/hAT-Tag1"
#12"DNA/hAT-Tip100"
#14"DNA/hAT-Tip100?"
#14"DNA?/hAT-Tip100?"
#15"DNA/Kolobok"
#16"DNA/Merlin"
#17"DNA/MULE-MuDR"
#18"DNA/PIF-Harbinger"
#19"DNA/PiggyBac"
#20"DNA/TcMar"
#21"DNA/TcMar-Mariner"
#22"DNA/TcMar-Pogo"
#23"DNA/TcMar-Tc1"
#24"DNA/TcMar-Tc2"
#25"DNA/TcMar-Tigger"
#26"RC/Helitron"

cat data/12_PDALseq/rm_annotation/DNA*.bed data/12_PDALseq/rm_annotation/DNA/* data/12_PDALseq/rm_annotation/DNA\?/* data/12_PDALseq/rm_annotation/RC/Helitron.bed  | bedtools sort > data/14_Repeat_masker_class_enrichment/01_compile_annotations/DNA_transposons.bed

#LTR
#"LTR"
#"LTR/ERV1"
#"LTR/ERVK"
#"LTR/ERVL"
#"LTR/ERVL-MaLR"
#"LTR/Gypsy"
cat data/12_PDALseq/rm_annotation/LTR.bed data/12_PDALseq/rm_annotation/LTR/* | bedtools sort > data/14_Repeat_masker_class_enrichment/01_compile_annotations/LTR.bed

#SINE
#"Retroposon/SVA"
#"SINE/5S-Deu-L2"
#"SINE/Alu"
#"SINE/MIR"
#"SINE/tRNA"
#"SINE/tRNA-Deu"
#"SINE/tRNA-RTE"
cat data/12_PDALseq/rm_annotation/Retroposon/SVA.bed data/12_PDALseq/rm_annotation/SINE/* | bedtools sort > data/14_Repeat_masker_class_enrichment/01_compile_annotations/SINE.bed

#LINE
#"LINE/CR1"
#"LINE/Dong-R4"
#"LINE/I-Jockey"
#"LINE/L1"
#"LINE/L1-Tx1"
#"LINE/L2"
#"LINE/Penelope"
#"LINE/RTE-BovB"
#"LINE/RTE-X"
cat  data/12_PDALseq/rm_annotation/LINE/* | bedtools sort > data/14_Repeat_masker_class_enrichment/01_compile_annotations/LINE.bed

#Psuedogene
#"rRNA"
#"scRNA"
#"snRNA"
#"srpRNA"
#"tRNA"
cat data/12_PDALseq/rm_annotation/*RNA* | bedtools sort > data/14_Repeat_masker_class_enrichment/01_compile_annotations/Psuedogene.bed

