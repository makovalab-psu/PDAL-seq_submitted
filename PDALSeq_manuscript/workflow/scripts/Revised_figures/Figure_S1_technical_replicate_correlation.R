library(tidyverse)
library(viridis)

df = read.csv("results/03_Correlation_and_black_list_unmappable_regions/Technical_replicate_correlation.csv")
df2 = read.csv("results/21_add_K562_to_analysis/Technical_replicate_correlation.csv")
df = bind_rows(df, df2 %>% select("Sample_1", "Sample_2", "Species",  "Window",   "Pearson",  "Spearman"))

df.plot = df[-which(grepl("TruSeq", df$Sample_1)),]
df.plot = df.plot[-which(grepl("TruSeq", df.plot$Sample_2)),]

unique(df.plot$Sample_1)

Sample_levels = c("Hsapien_GM24385_LCL_PDALseq_js1014A20_1_S2_R1_001",                                 
                  "Hsapien_GM24385_LCL_PDALseq_js1014A20_2_S3_R1_001",                                
                  "Hsapien_GM24385_LCL_PDALseq_js1014A40_1_S4_R1_001",                                 
                  "Hsapien_GM24385_LCL_PDALseq_js1014A40_2_S5_R1_001",                                 
                  "Hsapien_GM24143_LCL_PDALseq_js2007d20mM_S2_R1_001",                                 
                  "Hsapien_GM24143_LCL_PDALseq_js2007e20mM_S3_R1_001",                                 
                  "Hsapien_GM24143_LCL_PDALseq_js2007f40mM_S4_R1_001",                                 
                  "Hsapien_GM24143_LCL_PDALseq_js2007g40mM_S5_R1_001",                                 
                  "Hsapien_HEK293T_epithelial_PDALseq_js2008b20mM_S2_R1_001",                          
                  "Hsapien_HEK293T_epithelial_PDALseq_js2008c20mM_S3_R1_001",                          
                  "Hsapien_HEK293T_epithelial_PDALseq_js2008d40mM_S4_R1_001",                          
                  "Hsapien_HEK293T_epithelial_PDALseq_js2008e40mM_S5_R1_001",                          
                  "Hsapien_CCL86-Raji_LCL_PDALseq_ATCC_20mM_A_S2_R1_001",                              
                  "Hsapien_CCL86-Raji_LCL_PDALseq_ATCC_20mM_B_S3_R1_001",                              
                  "Hsapien_CCL86-Raji_LCL_PDALseq_ATCC_40mM_A_S4_R1_001",                              
                  "Hsapien_CCL86-Raji_LCL_PDALseq_ATCC_40mM_B_S5_R1_001",                              
                  "Hsapien_CCL86-Raji_LCL_PDALseq_F_1_20mM_A_S7_R1_001",                               
                  "Hsapien_CCL86-Raji_LCL_PDALseq_F_1_20mM_B_S8_R1_001",                               
                  "Hsapien_CCL86-Raji_LCL_PDALseq_F_1_40mM_A_S9_R1_001",                               
                  "Hsapien_CCL86-Raji_LCL_PDALseq_F_1_40mM_B_S10_R1_001",
                  "Hsa_K562_LCL_PDALSeq_js4003_20mM1_S2_R1_001",                                       
                  "Hsa_K562_LCL_PDALSeq_js4003_20mM2_S3_R1_001",                                       
                  "Hsa_K562_LCL_PDALSeq_js4003_40mM1_S4_R1_001",                                       
                  "Hsa_K562_LCL_PDALSeq_js4003_40mM2_S5_R1_001",                              
                  "Ptroglodytes_AG18354_LCL_PDALseq_js2016b20mM_S2_R1_001",                            
                  "Ptroglodytes_AG18354_LCL_PDALseq_js2016c20mM_S3_R1_001",                            
                  "Ptroglodytes_AG18354_LCL_PDALseq_js2016e40mM_S4_R1_001",                            
                  "Ptroglodytes_AG18354_LCL_PDALseq_js2029Mb20mM_S1_R1_001",                           
                  "Ptroglodytes_AG18354_LCL_PDALseq_js2029Mc20mM_S2_R1_001",                           
                  "Ptroglodytes_AG18354_LCL_PDALseq_js2029Md40mM_S3_R1_001",                           
                  "Ptroglodytes_AG18354_LCL_PDALseq_js2029Me40mM_S4_R1_001",                           
                  "Ptroglodytes_AG18358_LCL_PDALseq_js2029Fb20mM_S1_R1_001",                           
                  "Ptroglodytes_AG18358_LCL_PDALseq_js2029Fc20mM_S2_R1_001",                           
                  "Ptroglodytes_AG18358_LCL_PDALseq_js2029Fd40mM_S3_R1_001",                           
                  "Ptroglodytes_AG18358_LCL_PDALseq_js2029Fe40mM_S4_R1_001",                           
                  "Ppaniscus_KB8711_Fibroblast_PDALseq_js2025b20mM_S2_R1_001",                         
                  "Ppaniscus_KB8711_Fibroblast_PDALseq_js2025c20mM_S3_R1_001",                         
                  "Ppaniscus_KB8711_Fibroblast_PDALseq_js2025d40mM_S4_R1_001",                         
                  "Ppaniscus_KB8711_Fibroblast_PDALseq_js2025e40mM_S5_R1_001",                         
                  "Ppaniscus_KB14571_Fibroblast_PDALseq_PPA_20_1_1_241205_S1_R1_001",                  
                  "Ppaniscus_KB14571_Fibroblast_PDALseq_PPA_20_2_2_241205_S2_R1_001",                  
                  "Ppaniscus_KB14571_Fibroblast_PDALseq_PPA_40_1_3_241205_S3_R1_001",                  
                  "Ppaniscus_KB14571_Fibroblast_PDALseq_PPA_40_1_test_5_241205_S5_R1_001",             
                  "Ppaniscus_KB14571_Fibroblast_PDALseq_PPA_40_2_4_241205_S4_R1_001",                  
                  "Ggorilla_KB3781_Fibroblast_PDALseq_js2019c20mM_S2_R1_001",                          
                  "Ggorilla_KB3781_Fibroblast_PDALseq_js2022c20mM_S3_R1_001",                          
                  "Ggorilla_KB3781_Fibroblast_PDALseq_js2022d40mM_S4_R1_001",                          
                  "Ggorilla_KB3781_Fibroblast_PDALseq_js2022e40mM_S5_R1_001",                          
                  "Ggorilla_F7549_Fibroblast_PDALseq_GGO_20_1_1_241206_S1_R1_001",                     
                  "Ggorilla_F7549_Fibroblast_PDALseq_GGO_20_2_2_241206_S2_R1_001",                     
                  "Ggorilla_F7549_Fibroblast_PDALseq_GGO_20_2_cntl_5_241206_S5_R1_001",                
                  "Ggorilla_F7549_Fibroblast_PDALseq_GGO_40_1_3_241206_S3_R1_001",                     
                  "Ggorilla_F7549_Fibroblast_PDALseq_GGO_40_2_4_241206_S4_R1_001",                     
                  "Pabelii_AG06213_Fibroblast_PDALseq_js3008b_20241211_SOR-AG06213-20mM-1-2_S2_R1_001",
                  "Pabelii_AG06213_Fibroblast_PDALseq_js3008c_20241211_SOR-AG06213-20mM-2-3_S3_R1_001",
                  "Pabelii_AG06213_Fibroblast_PDALseq_js3008d_20241211_SOR-AG06213-40mM-1-4_S4_R1_001",
                  "Pabelii_AG06213_Fibroblast_PDALseq_js3008e_20241211_SOR-AG06213-40mM-2-5_S5_R1_001",
                  "Ppygmaeus_AG05252_Fibroblast_PDALseq_20_1_240822_S1_R1_001",                        
                  "Ppygmaeus_AG05252_Fibroblast_PDALseq_20_2_240822_S2_R1_001",                        
                  "Ppygmaeus_AG05252_Fibroblast_PDALseq_40_1_240822_S3_R1_001",                        
                  "Ppygmaeus_AG05252_Fibroblast_PDALseq_40_2_240822_S4_R1_001",                        
                  "Ssyndactylus_JAMBI_LCL_PDALseq_js2012h20mM_S3_R1_001",                              
                  "Ssyndactylus_JAMBI_LCL_PDALseq_js2012j40mM_S4_R1_001",                              
                  "Ssyndactylus_JAMBI_LCL_PDALseq_js2016g20mM_S2_R1_001",                              
                  "Ssyndactylus_JAMBI_LCL_PDALseq_js2034c20mM_S1_R1_001",                              
                  "Ssyndactylus_JAMBI_LCL_PDALseq_js2034d40mM_S2_R1_001",                              
                  "Ssyndactylus_JAMBI_LCL_PDALseq_js2034e40mM_S3_R1_001",                              
                  "Ssyndactylus_Karenina_LCL_PDALseq_js3008g_20241218_SSY-Ka-20mM-1-2_S2_R1_001",      
                  "Ssyndactylus_Karenina_LCL_PDALseq_js3008h_20241218_SSY-Ka-20mM-2-3_S3_R1_001",      
                  "Ssyndactylus_Karenina_LCL_PDALseq_js3008i_20241218_SSY-Ka-40mM-1-4_S4_R1_001",      
                  "Ssyndactylus_Karenina_LCL_PDALseq_js3008j_20241218_SSY-Ka-40mM-2-5_S5_R1_001")

df.plot$Sample_1 = factor(df.plot$Sample_1, levels = Sample_levels)
df.plot$Sample_2 = factor(df.plot$Sample_2, levels = Sample_levels)


P = ggplot(df.plot, aes(x = Sample_1, y = Sample_2, fill = Spearman)) +
	facet_wrap(~Window) +	
	geom_tile() +
	theme_classic() +
	theme(axis.text.x = element_text(angle =90),
	      axis.text.y = element_blank(),
	      axis.title = element_blank()) +
	scale_fill_viridis(limits = c(0,1), option = "turbo") +
	coord_fixed()
	
ggsave("revised_plots/Figure_S1_technical_replicate_correlation.svg", P, width = 5, height = 10, scale = 3)

#fastq-dump --accession SRR11032812 --split-files --split-spot --gzip
