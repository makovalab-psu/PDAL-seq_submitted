#Snakefile for the siamang chr12/chr19 rerun (step 32): the main Snakefile plus step 17,
#which is commented out there. Run with -s workflow/siamang_fix.smk and explicit targets
#(see workflow/scripts/32_siamang_chr12_chr19_fix/README.md)

include: "Snakefile"
include: "rules/17_PDALseq_reads_in_CenSat.smk"
