
#write bed files for calculating enrichment

#SST1
grep SST1 resources/repeats/linneas_annotations/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SST1/Hsapien.bed
grep SST1 resources/repeats/linneas_annotations/mGorGor1.pri.cur.20231122.RepeatMasker-combo3.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SST1/Ggorilla.bed
grep SST1 resources/repeats/linneas_annotations/mPanPan1.pri.cur.20231122.RepeatMasker-combo3.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SST1/Ppaniscus.bed
grep SST1 resources/repeats/linneas_annotations/mPanTro3.pri.cur.20231122.RepeatMasker-combo.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SST1/Ptroglodytes.bed
grep SST1 resources/repeats/linneas_annotations/mPonAbe1.pri.cur.20231205.RepeatMasker-combo.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SST1/Pabelii.bed
grep SST1 resources/repeats/linneas_annotations/mPonPyg2.pri.cur.20231122.RepeatMasker-combo.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SST1/Ppygmaeus.bed
grep SST1 resources/repeats/linneas_annotations/mSymSyn1.pri.cur.20231205.RepeatMasker-combo.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | awk -v pair=chr12,chr19 -f workflow/scripts/swap_chrom_pair.awk | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SST1/Ssyndactylus.bed

#SATR1
grep SATR1 resources/repeats/linneas_annotations/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SATR1/Hsapien.bed
grep SATR1 resources/repeats/linneas_annotations/mGorGor1.pri.cur.20231122.RepeatMasker-combo3.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SATR1/Ggorilla.bed
grep SATR1 resources/repeats/linneas_annotations/mPanPan1.pri.cur.20231122.RepeatMasker-combo3.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SATR1/Ppaniscus.bed
grep SATR1 resources/repeats/linneas_annotations/mPanTro3.pri.cur.20231122.RepeatMasker-combo.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SATR1/Ptroglodytes.bed
grep SATR1 resources/repeats/linneas_annotations/mPonAbe1.pri.cur.20231205.RepeatMasker-combo.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SATR1/Pabelii.bed
grep SATR1 resources/repeats/linneas_annotations/mPonPyg2.pri.cur.20231122.RepeatMasker-combo.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SATR1/Ppygmaeus.bed
grep SATR1 resources/repeats/linneas_annotations/mSymSyn1.pri.cur.20231205.RepeatMasker-combo.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | awk -v pair=chr12,chr19 -f workflow/scripts/swap_chrom_pair.awk | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SATR1/Ssyndactylus.bed

#SATR2
grep SATR2 resources/repeats/linneas_annotations/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SATR2/Hsapien.bed
grep SATR2 resources/repeats/linneas_annotations/mGorGor1.pri.cur.20231122.RepeatMasker-combo3.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SATR2/Ggorilla.bed
grep SATR2 resources/repeats/linneas_annotations/mPanPan1.pri.cur.20231122.RepeatMasker-combo3.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SATR2/Ppaniscus.bed
grep SATR2 resources/repeats/linneas_annotations/mPanTro3.pri.cur.20231122.RepeatMasker-combo.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SATR2/Ptroglodytes.bed
grep SATR2 resources/repeats/linneas_annotations/mPonAbe1.pri.cur.20231205.RepeatMasker-combo.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SATR2/Pabelii.bed
grep SATR2 resources/repeats/linneas_annotations/mPonPyg2.pri.cur.20231122.RepeatMasker-combo.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SATR2/Ppygmaeus.bed
grep SATR2 resources/repeats/linneas_annotations/mSymSyn1.pri.cur.20231205.RepeatMasker-combo.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | awk -v pair=chr12,chr19 -f workflow/scripts/swap_chrom_pair.awk | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/SATR2/Ssyndactylus.bed

#LSAU
grep LSAU resources/repeats/linneas_annotations/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/LSAU/Hsapien.bed
grep LSAU resources/repeats/linneas_annotations/mGorGor1.pri.cur.20231122.RepeatMasker-combo3.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/LSAU/Ggorilla.bed
grep LSAU resources/repeats/linneas_annotations/mPanPan1.pri.cur.20231122.RepeatMasker-combo3.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/LSAU/Ppaniscus.bed
grep LSAU resources/repeats/linneas_annotations/mPanTro3.pri.cur.20231122.RepeatMasker-combo.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/LSAU/Ptroglodytes.bed
grep LSAU resources/repeats/linneas_annotations/mPonAbe1.pri.cur.20231205.RepeatMasker-combo.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/LSAU/Pabelii.bed
grep LSAU resources/repeats/linneas_annotations/mPonPyg2.pri.cur.20231122.RepeatMasker-combo.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/LSAU/Ppygmaeus.bed
grep LSAU resources/repeats/linneas_annotations/mSymSyn1.pri.cur.20231205.RepeatMasker-combo.out | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | awk -v pair=chr12,chr19 -f workflow/scripts/swap_chrom_pair.awk | bedtools sort | bedtools merge > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/LSAU/Ssyndactylus.bed


