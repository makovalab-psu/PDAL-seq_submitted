

#echo $(cat $1 | awk '{print $1 "\t" $2 "\t" $3}' > $2)

declare -a nonb=("APR" "CRU" "DR" "GQ" "IR" "MR" "SLS" "STR" "TRI" "Z")

## now loop through the above array
for i in "${nonb[@]}"
do
   echo "$i"
   #Hsapien
   cat resources/nonb_annotations/chm13v2.0.$i.bed | awk '{print $1 "\t" $2 "\t" $3}' > results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/Hsapien/$i.bed
   #Ptroglodytes
   cat resources/nonb_annotations/mPanTro3.pri.cur.20231122.$i.bed | awk '{print $1 "\t" $2 "\t" $3}' > results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/Ptroglodytes/$i.bed
   #Ppaniscus
   cat resources/nonb_annotations/mPanPan1.pri.cur.20231122.$i.bed | awk '{print $1 "\t" $2 "\t" $3}' > results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/Ppaniscus/$i.bed
   #Ggorilla
   cat resources/nonb_annotations/mGorGor1.pri.cur.20231122.$i.bed | awk '{print $1 "\t" $2 "\t" $3}' > results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/Ggorilla/$i.bed
   #Pabelii
   cat resources/nonb_annotations/mPonAbe1.pri.cur.20231205.$i.bed | awk '{print $1 "\t" $2 "\t" $3}' > results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/Pabelii/$i.bed
   #Ppygmaeus
   cat resources/nonb_annotations/mPonPyg2.pri.cur.20231122.$i.bed | awk '{print $1 "\t" $2 "\t" $3}' > results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/Ppygmaeus/$i.bed
   #Ssyndactylus (annotation is on mSymSyn1 v2.1; swap chr12/chr19 back to the v2.0 labels of the reference)
   cat resources/nonb_annotations/mSymSyn1.pri.cur.20240514.$i.bed | awk '{print $1 "\t" $2 "\t" $3}' | awk -v pair=chr12,chr19 -f workflow/scripts/swap_chrom_pair.awk > results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/Ssyndactylus/$i.bed
done

echo "done"

