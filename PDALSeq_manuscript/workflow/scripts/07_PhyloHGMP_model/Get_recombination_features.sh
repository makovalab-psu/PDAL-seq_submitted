#Bed format the different types of recombination features
grep -v "#" resources/functional_annotations/Refseq/genomic_ucsc_names.gff | awk '$3 == "meiotic_recombination_region" {print $1 "\t" $4 "\t" $5 "\t" $3}' > resources/functional_annotations/Refseq/meiotic_recombination_region.bed
grep -v "#" resources/functional_annotations/Refseq/genomic_ucsc_names.gff | awk '$3 == "mitotic_recombination_region" {print $1 "\t" $4 "\t" $5 "\t" $3}' > resources/functional_annotations/Refseq/mitotic_recombination_region.bed
grep -v "#" resources/functional_annotations/Refseq/genomic_ucsc_names.gff | awk '$3 == "non_allelic_homologous_recombination_region" {print $1 "\t" $4 "\t" $5 "\t" $3}' > resources/functional_annotations/Refseq/non_allelic_homologous_recombination_region.bed
grep -v "#" resources/functional_annotations/Refseq/genomic_ucsc_names.gff | awk '$3 == "chromosome_breakpoint" {print $1 "\t" $4 "\t" $5 "\t" $3}' > resources/functional_annotations/Refseq/chromosome_breakpoint.bed
grep -v "#" resources/functional_annotations/Refseq/genomic_ucsc_names.gff | awk '$3 == "recombination_feature" {print $1 "\t" $4 "\t" $5 "\t" $3}' > resources/functional_annotations/Refseq/recombination_feature.bed

#Combine annotations

cat resources/functional_annotations/Refseq/meiotic_recombination_region.bed resources/functional_annotations/Refseq/mitotic_recombination_region.bed resources/functional_annotations/Refseq/non_allelic_homologous_recombination_region.bed resources/functional_annotations/Refseq/chromosome_breakpoint.bed resources/functional_annotations/Refseq/recombination_feature.bed | bedtools sort | bedtools merge > resources/functional_annotations/Saswats_annotations/Hsapien/Recombination.bed 
