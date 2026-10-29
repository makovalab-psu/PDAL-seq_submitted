# Figures and datasets

Which script makes each figure and SI dataset in the published paper. Scripts are in
`workflow/scripts/` and are run from `PDALSeq_manuscript/` after the workflow (see
[README](README.md#figures-and-tables)).

Numbering follows the final manuscript and SI Appendix. Most script and output file names use
the numbering of earlier drafts, so **the number in a script's name is usually not the figure
it makes**. For example, `Figure_S6_general_linear_model…R` makes Fig. S6, but
`Figure_S14_PCA…R` makes Fig. S9 and `new_plots/Figure_6_Sat_arrays.R` makes Fig. 5.
Scripts not listed here are earlier drafts that are not used in the paper.

## Main figures

| Figure | Content | Script | Data from steps |
|---|---|---|---|
| 1 | Experimental design | none (drawn by hand) | – |
| 2A–C | IGV traces at the MYC, KRAS and ZMIZ2 promoters (HEK293T) | none (IGV screenshots) | 01, 08 |
| 2D–E | Correlation with non-B motif coverage; reads in syntenic vs. non-syntenic regions | `Revised_figures/Figure_2_Genome_wide_comparison_to_nonb_motifs_revised.R` | 04, 19, 21 |
| 3 | Ape HMM states (A), human chromosome ideogram of states (B), signal around core/stochastic origins (C) | `Revised_figures/Figure_3_Ape_PDAL-Seq_MVGHMM.R` | 07, 18, 21 |
| 4 | Human HMM states (A) and their repeat content (B) | `Revised_figures/Figure_4_Human_HMM.R` | 12, 26, 27 |
| 5A–C | Satellite enrichment and non-B coverage (A), HSat3 (AATGG)n reference-free counts (B), CD spectra (C) | `new_plots/Figure_6_Sat_arrays.R` | 06, 17; CD data not from the workflow |
| 5D | G-zipper hairpin model | none (drawn by hand) | – |

## SI figures

| Figure | Content | Script | Data from steps |
|---|---|---|---|
| S1 | Technical-replicate correlation | `Revised_figures/Figure_S1_technical_replicate_correlation.R` | 03, 21, 25 |
| S2 | HOXD11 promoter | none (IGV screenshot) | 01, 08 |
| S3 | Genome-wide PDAL-Seq reads, non-B motif and HSat2/3 density per species | `Revised_figures/Figure_S2_Genome_wide_plots.R` | 04 |
| S4 | Reads vs. non-B motif coverage in 1-Mbp windows, by motif | `Revised_figures/Figure_S3_Genome_wide_correlation_plots.R` | 04, 21 |
| S5 | Strand-specific coverage around non-B motifs | `Revised_figures/Figure_S4_strand_specific_analysis.R` | 20 |
| S6 | Poisson GLM parameters, unique 1-kbp windows | `Revised_figures/Figure_S6_general_linear_model_analysis_of_reads_in_unique_1kb_windows.R` | 16, 21 |
| S7 | Comparison with G4/R-loop ChIP-seq and CUT&Tag (HEK293T) | `Figures_Tables/Figure_S5_CHIP_CutRun_Benchmarking.R` | 08 |
| S8 | Transition probabilities of the 8-state ape HMM | `Figures_Tables/Figure_S11_HMM_transition_probabilities.R` | 07 |
| S9 | PCA of PDAL-Seq and TruSeq in Hominoidae alignment blocks | `Revised_figures/Figure_S14_PCA_15_cell_lines_7_way_alignment.R` | 26 |
| S10 | Lymphoblast (A) and fibroblast (B) HMMs | `29_Lymphoblast_specific_HMM/Figure_SI_Lymphoblast_HMM.R` and `28_Fibroblast_specific_HMM/Figure_SI_Fibroblast_HMM.R` (run by steps 29 and 28) | 28, 29 |
| S11 | CpG 5mC levels by ape HMM state | `new_plots/Figure_4_Methylation_origins_TSS_in_states.R` (the 5mC heatmap, `P.5mC`) | 10 |
| S12 | Signal around alternative replication origins | `Revised_figures/Figure_S12_PDALSeq_signal_around_alternative_origin_sets.R` | 22 |
| S13 | Signal around TSSs | `Revised_figures/Figure_S13_TSS_plot.R` | 12, 21 |
| S14 | PDAL-Seq and TruSeq enrichment in ape satellite arrays | `new_plots/Figure_S17_TruSeq_Sat_coverage.R` | 17 |
| S15 | Siamang subterminal alpha satellites | `Figures_Tables/Figure_S4_Analysis_of_PDALseq_reads_in_subterminal_Siamang_alpha_sats.R` | 04 |
| S16 | PDAL-Seq, DR motifs and HSat2/3 in homologous great ape chromosomes | `Figures_Tables/Figure_SX_Hsat_ideograms.R` | 04 |
| S17A–B | PCA of CD spectra (A); SEC-UV of (AATGG)n (B) | `new_plots/Figure_6_Sat_arrays.R` (`PE` and `P.chrom`) | not from the workflow (biophysical data) |
| S17C | Native PAGE gel | none (gel image) | – |
| S17D | SEC-MALS molecular weights | none (see Dataset S11) | – |
| S18 | Cell-cycle (DMSO arrest) contribution: replicate correlation (B), fold change (C), window classes (D); panel A is flow cytometry | `Revised_figures/Figure_S7_S-phase_contribution_to_PDAL-Seq_signal.R` | 25 |
| S19 | Cell-cycle-specific regions vs. non-B motifs, functional elements and ape HMM states | `Revised_figures/Figure_S8_analysis_on_functional_element_contribution_to_cell_cycle_effects.R` | 25 |
| S20 | PCA of alignment blocks for each species group | `Figures_Tables/Figure_S8_PDALseq_PCA_different_taxon.R` | 07 |
| S21 | Correlation heatmap of signal in Hominoidae alignment blocks | `new_plots/Figure_SX_Correllation_clustering_on_PDAL-Seq_reads_in_allignment.R` | 07 |
| S22 | Pairwise correlations by comparison type (Kruskal-Wallis) | `Figures_Tables/Figure_SX_correllation_between_male_female_transformed_non_transformed.R` | 07 |
| S23 | K-means (A) and HMM BIC (B) scree plots for each species group | `Figures_Tables/Figure_S9_PDALseq_screeplot_different_taxon.R` | 07 |
| S24 | Ape HMM state means at 2, 8, 18 and 30 states, and BIC | `Figures_Tables/Figure_S10_BIC_versus_hidden_states_interpratation.R` | 07 |
| S25 | Scree plots for the cell-type HMMs (A) and the human HMM (B) | `Revised_figures/Figure_S24_Scree_plots_for_HMMs.R` (its `ggsave` call has no file name; save the plot `P` by hand) | 26, 28, 29 |
| S26 | S1 digestion test gel | none (gel image) | – |

Scripts in `Figures_Tables/` that take the output path as an argument are run by the rules in
`workflow/rules/Figures_Tables.smk`, e.g. `Rscript workflow/scripts/Figures_Tables/Figure_S11_HMM_transition_probabilities.R out.svg`.

## SI datasets

| Dataset | Content | Script | Data from steps |
|---|---|---|---|
| S1 | Cell lines | none (written by hand) | – |
| S2 | Mapping statistics per technical replicate | `Figures_Tables/Table_S1_preprocessing_statistics.R`; K562 and DMSO rows from `30_update_SI_tables_with_revisions_data/Compile_Table_S1.R` | 02, 30 |
| S3 | Mapping statistics per merged cell line | `Figures_Tables/Table_S2_preprocessing_statistics_by_sample.R`; K562 and DMSO rows from `30_…/Compile_Table_S2.R` | 02, 30 |
| S4 | Multi-mapping regions and unique mapping rates | `Figures_Tables/Table_S3_unique_maps_inside_and_outside_blacklist.R`; K562 and DMSO rows from `30_…/Compile_Table_S3.R` | 03, 30 |
| S5 | McFadden's R² of Poisson GLMs, unique 1-kbp windows | `16_lm_analysis_of_nonbmotifs/glm_results_1kbp_unique_mapping_windows.R` (writes `Global_fit_R2.csv`; run by hand) | 16 |
| S6 | GLM coefficients and VIFs, unique 1-kbp windows | same script (writes `Global_fit_coefficients.csv`) | 16 |
| S7 | GO and KEGG enrichment of human HMM state genes | `Figures_Tables/Figure_4_Analysis_of_PDAL-seq_reads.R` (writes `tables/Table_SX_Human_state_8_Enrichr_analysis.csv`) | 09 |
| S8 | Homologous chromosomes kept by synteny filtering | `Process_alignments/workflow/scripts/identify_synteny.R` (output in `Process_alignments/data/11_identify_syteny/`) | Process_alignments |
| S9 | Alignment filtering statistics | `Process_alignments/workflow/scripts/Compile_alignment_stats.R` and `Median_allinment_length.R` | Process_alignments |
| S10 | Unique mapping rates inside alignment blocks | `Figures_Tables/Table_S7_unique_mapping_rates_in_alingments.R` | 07 |
| S11 | SEC-MALS molecular weights of (AATGG)n | none (instrument software) | – |
| S12 | Pairwise correlations in Hominoidae alignment blocks | `Figures_Tables/Table_SX_PDALseq_correlation_in_alignments.R` | 07 |

Spreadsheet layout and column descriptions in the published datasets were added by hand.
