####################################################################
# Step 31 — UCSC track hubs
#
# Builds one UCSC track hub per ape T2T genome (7 hubs). Each hub has:
#   - PDAL-Seq signal for every PDAL-Seq cell line of that species
#     (bigWig, CPM-normalized coverage from the merged All_reads BAMs;
#     human also includes K562 from step 21)
#   - Ape HMM states (bigBed 9): the 8-state hominoidae multivariate
#     Gaussian HMM from step 07, projected onto each species through the
#     shared row order of the alignment BEDs, with the state labels and
#     colors of Fig. 3 (Revised_figures/Figure_3_Ape_PDAL-Seq_MVGHMM.R)
#   - Non-B DNA motif annotations used in the manuscript (bigBed 6, one
#     track per motif type; G4 keeps its strand)
#   - Human only: published HEK293T G4 CUT&Tag data mapped in step 08
#     (bigWig, CPM)
#
# Human uses UCSC hs1; the other apes use their UCSC GenArk assemblies.
# Chromosome names are converted to each assembly's native names with
# the UCSC chromAlias file. Records on sequences missing from chromAlias
# are dropped and counted in the log.
#
# Upstream BAMs, state calls and alignment BEDs are ancient() so steps
# 01, 07, 08 and 21 are never re-triggered.
# Output: results/31_UCSC_track_hubs/{species}/  (hub.txt + data files).
# Copy each folder to a public web server and load its hub.txt URL in
# the UCSC Genome Browser (My Data > Track Hubs).
####################################################################

HUB_EMAIL_31 = "kdm16@psu.edu"

_UCSC_31 = "https://hgdownload.soe.ucsc.edu"

# species: (UCSC genome, common name, assembly name, non-B file prefix, Fig. 2D species color)
_SPECIES_31 = {
    "Hsapien":      ("hs1",             "human",              "CHM13v2.0", "chm13v2.0",                 "158,158,158"),
    "Ptroglodytes": ("GCA_028858775.2", "chimpanzee",         "mPanTro3",  "mPanTro3.pri.cur.20231122", "190,36,47"),
    "Ppaniscus":    ("GCA_029289425.2", "bonobo",             "mPanPan1",  "mPanPan1.pri.cur.20231122", "220,147,182"),
    "Ggorilla":     ("GCA_029281585.2", "gorilla",            "mGorGor1",  "mGorGor1.pri.cur.20231122", "28,125,144"),
    "Ppygmaeus":    ("GCA_028885625.2", "Bornean orangutan",  "mPonPyg2",  "mPonPyg2.pri.cur.20231122", "242,217,131"),
    "Pabelii":      ("GCA_028885655.2", "Sumatran orangutan", "mPonAbe1",  "mPonAbe1.pri.cur.20231205", "250,150,50"),
    "Ssyndactylus": ("GCA_028878055.2", "siamang",            "mSymSyn1",  "mSymSyn1.pri.cur.20240514", "49,55,87"),
}

# Chromosome pairs to swap in the assembly-name column (2) of a UCSC chromAlias file.
# GCA_028878055.2 names the 98 Mb siamang chromosome (CM054520.2) chr12_hap1 and the 146 Mb one
# (CM054527.2) chr19_hap1; mSymSyn1 v2.1, used for our FASTA, BAMs, alignments and non-B calls,
# swapped those two labels.
_ALIAS_CHROM_SWAP_31 = {
    "GCA_028878055.2": "chr12,chr19",
}

# non-B type in file names: (track label, description)
_NONB_31 = {
    "APR": ("APR", "A-phased repeats"),
    "CRU": ("CRU", "cruciform motifs"),
    "DR":  ("DR",  "direct repeats"),
    "GQ":  ("G4",  "G-quadruplex motifs (Quadron)"),
    "IR":  ("IR",  "inverted repeats"),
    "MR":  ("MR",  "mirror repeats"),
    "SLS": ("SLS", "slipped-strand motifs"),
    "STR": ("STR", "short tandem repeats"),
    "TRI": ("TRI", "triplex motifs"),
    "Z":   ("Z",   "Z-DNA motifs"),
}

_K562_BAM_31 = "data/21_add_K562_to_analysis/01_preprocess/All_reads/Hsapien_K562_PDALseq.bam"

# PDAL-Seq cell lines per species, e.g. {"Hsapien": ["GM24385", ...]}
_PDALSEQ_31 = {s: [] for s in _SPECIES_31}
for _exp in sorted(DATA):
    _species, _biosample, _technique = _exp.split("_")
    if _technique == "PDALseq":
        _PDALSEQ_31[_species].append(_biosample)
_PDALSEQ_31["Hsapien"].append("K562")

_G4_CUTTAG_31 = df_benchmark[df_benchmark["Assay"] == "G4 CUT&Tag"]

_HUB_SPECIES_RE_31 = "|".join(_SPECIES_31)


def _genome_31(wildcards):
    return _SPECIES_31[wildcards.hub_species][0]


def _ucsc_url_31(genome, suffix):
    if genome == "hs1":
        return f"{_UCSC_31}/goldenPath/hs1/bigZips/hs1.{suffix}"
    acc = genome.replace("_", "")
    path = "/".join([genome[:3], acc[3:6], acc[6:9], acc[9:12], genome])
    return f"{_UCSC_31}/hubs/{path}/{genome}.{suffix}"


def _pdalseq_bam_31(wildcards):
    if wildcards.biosample == "K562":
        return _K562_BAM_31
    return f"data/01_preprocess_data/All_reads/{wildcards.hub_species}_{wildcards.biosample}_PDALseq.bam"


def _hub_files_31(wildcards):
    sp = wildcards.hub_species
    files = [f"results/31_UCSC_track_hubs/{sp}/pdalseq_{b}.bw" for b in _PDALSEQ_31[sp]]
    files.append(f"results/31_UCSC_track_hubs/{sp}/ape_hmm_states.bb")
    files += [f"results/31_UCSC_track_hubs/{sp}/nonb_{n}.bb" for n in _NONB_31]
    if sp == "Hsapien":
        files += [f"results/31_UCSC_track_hubs/{sp}/g4_cuttag_{r}.bw" for r in _G4_CUTTAG_31["Run"]]
    return files


rule download_ucsc_genome_files_31:
    output:
        alias = "data/31_UCSC_track_hubs/00_ucsc_genome_files/{genome}.chromAlias.txt",
        sizes = "data/31_UCSC_track_hubs/00_ucsc_genome_files/{genome}.chrom.sizes"
    params:
        alias_url = lambda wc: _ucsc_url_31(wc.genome, "chromAlias.txt"),
        sizes_url = lambda wc: _ucsc_url_31(wc.genome, "chrom.sizes.txt"),
        swap = lambda wc: _ALIAS_CHROM_SWAP_31.get(wc.genome, "")
    log:
        "logs/31_UCSC_track_hubs/00_ucsc_genome_files/{genome}.txt"
    shell:
        "curl -fsSL {params.alias_url} 2> {log} "
        "| awk -v pair='{params.swap}' -v col=2 -f workflow/scripts/swap_chrom_pair.awk > {output.alias} && "
        "curl -fsSL -o {output.sizes} {params.sizes_url} 2>> {log}"


rule pdalseq_bigwig_31:
    wildcard_constraints:
        hub_species = _HUB_SPECIES_RE_31,
        biosample = r"[A-Za-z0-9-]+"
    input:
        bam = lambda wc: ancient(_pdalseq_bam_31(wc)),
        alias = lambda wc: f"data/31_UCSC_track_hubs/00_ucsc_genome_files/{_genome_31(wc)}.chromAlias.txt",
        sizes = lambda wc: f"data/31_UCSC_track_hubs/00_ucsc_genome_files/{_genome_31(wc)}.chrom.sizes"
    output:
        "results/31_UCSC_track_hubs/{hub_species}/pdalseq_{biosample}.bw"
    log:
        "logs/31_UCSC_track_hubs/01_pdalseq_bigwig/{hub_species}_{biosample}.txt"
    conda:
        "../env/31_track_hubs.yml"
    resources:
        mem_mb=12000
    shell:
        "bash workflow/scripts/31_UCSC_track_hubs/bam_to_cpm_bigwig.sh "
        "{input.bam} {input.alias} {input.sizes} {output} 2> {log}"


rule g4_cuttag_bigwig_31:
    input:
        bam = ancient("data/08_CHIP_CutRun_Benchmarking/04_remove_dup/{sra}.bam"),
        alias = "data/31_UCSC_track_hubs/00_ucsc_genome_files/hs1.chromAlias.txt",
        sizes = "data/31_UCSC_track_hubs/00_ucsc_genome_files/hs1.chrom.sizes"
    output:
        "results/31_UCSC_track_hubs/Hsapien/g4_cuttag_{sra}.bw"
    log:
        "logs/31_UCSC_track_hubs/02_g4_cuttag_bigwig/{sra}.txt"
    conda:
        "../env/31_track_hubs.yml"
    resources:
        mem_mb=12000
    shell:
        "bash workflow/scripts/31_UCSC_track_hubs/bam_to_cpm_bigwig.sh "
        "{input.bam} {input.alias} {input.sizes} {output} 2> {log}"


rule ape_hmm_bigbed_31:
    wildcard_constraints:
        hub_species = _HUB_SPECIES_RE_31
    input:
        states = ancient("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/GHMM_states.csv"),
        blocks = lambda wc: ancient(f"resources/alignments/hominoidae/{wc.hub_species}.bed"),
        alias = lambda wc: f"data/31_UCSC_track_hubs/00_ucsc_genome_files/{_genome_31(wc)}.chromAlias.txt",
        sizes = lambda wc: f"data/31_UCSC_track_hubs/00_ucsc_genome_files/{_genome_31(wc)}.chrom.sizes"
    output:
        bed = temp("data/31_UCSC_track_hubs/03_ape_hmm_states/{hub_species}.bed"),
        bb = "results/31_UCSC_track_hubs/{hub_species}/ape_hmm_states.bb"
    log:
        "logs/31_UCSC_track_hubs/03_ape_hmm_states/{hub_species}.txt"
    conda:
        "../env/31_track_hubs.yml"
    resources:
        mem_mb=8000
    shell:
        "(python workflow/scripts/31_UCSC_track_hubs/hmm_states_to_bed9.py {input.states} {input.blocks} "
        "| python workflow/scripts/31_UCSC_track_hubs/rename_chroms.py {input.alias} "
        "| LC_ALL=C sort -k1,1 -k2,2n > {output.bed}) 2> {log} && "
        "bedToBigBed -type=bed9 {output.bed} {input.sizes} {output.bb} 2>> {log}"


rule nonb_bigbed_31:
    wildcard_constraints:
        hub_species = _HUB_SPECIES_RE_31,
        nonb = "|".join(_NONB_31)
    input:
        nonb = lambda wc: f"resources/nonb_annotations/{_SPECIES_31[wc.hub_species][3]}.{wc.nonb}.bed",
        alias = lambda wc: f"data/31_UCSC_track_hubs/00_ucsc_genome_files/{_genome_31(wc)}.chromAlias.txt",
        sizes = lambda wc: f"data/31_UCSC_track_hubs/00_ucsc_genome_files/{_genome_31(wc)}.chrom.sizes"
    output:
        bed = temp("data/31_UCSC_track_hubs/04_nonb/{hub_species}/{nonb}.bed"),
        bb = "results/31_UCSC_track_hubs/{hub_species}/nonb_{nonb}.bb"
    params:
        label = lambda wc: _NONB_31[wc.nonb][0]
    log:
        "logs/31_UCSC_track_hubs/04_nonb/{hub_species}/{nonb}.txt"
    conda:
        "../env/31_track_hubs.yml"
    resources:
        mem_mb=8000
    shell:
        "(bash workflow/scripts/31_UCSC_track_hubs/nonb_to_bed6.sh {input.nonb} {params.label} "
        "| python workflow/scripts/31_UCSC_track_hubs/rename_chroms.py {input.alias} "
        "| LC_ALL=C sort -k1,1 -k2,2n > {output.bed}) 2> {log} && "
        "bedToBigBed -type=bed6 {output.bed} {input.sizes} {output.bb} 2>> {log}"


rule hub_txt_31:
    wildcard_constraints:
        hub_species = _HUB_SPECIES_RE_31
    input:
        _hub_files_31
    output:
        "results/31_UCSC_track_hubs/{hub_species}/hub.txt"
    run:
        sp = wildcards.hub_species
        genome, common, assembly, _, color = _SPECIES_31[sp]

        def stanza(lines, indent=""):
            return "".join(f"{indent}{line}\n" for line in lines) + "\n"

        out = stanza([
            f"hub PDALSeq_{sp}",
            f"shortLabel PDAL-Seq {genome if genome == 'hs1' else assembly}",
            f"longLabel PDAL-Seq ssDNA, ape HMM states, non-B DNA motifs: {common} ({assembly})",
            "useOneFile on",
            f"email {HUB_EMAIL_31}",
        ])
        out += stanza([f"genome {genome}"])

        out += stanza([
            "track pdalseq",
            "compositeTrack on",
            "shortLabel PDAL-Seq",
            "longLabel PDAL-Seq signal (CPM, all deduplicated reads)",
            "type bigWig",
            "visibility full",
            "autoScale on",
            "maxHeightPixels 100:40:8",
        ])
        for b in _PDALSEQ_31[sp]:
            out += stanza([
                f"track pdalseq_{b.replace('-', '_')}",
                "parent pdalseq on",
                f"shortLabel {b[:17]}",
                f"longLabel PDAL-Seq, {common} {b} (CPM)",
                "type bigWig",
                f"bigDataUrl pdalseq_{b}.bw",
                f"color {color}",
            ], indent="    ")

        out += stanza([
            "track ape_hmm_states",
            "shortLabel Ape HMM states",
            "longLabel PDAL-Seq ape HMM states (8-state Gaussian HMM, Hominoidae alignment blocks)",
            "type bigBed 9",
            "itemRgb on",
            "visibility dense",
            "bigDataUrl ape_hmm_states.bb",
        ])

        out += stanza([
            "track nonb",
            "compositeTrack on",
            "shortLabel Non-B DNA motifs",
            "longLabel Predicted non-B DNA motifs used in the PDAL-Seq manuscript (GFA; G4 by Quadron)",
            "type bigBed 6",
            "visibility dense",
        ])
        for n, (label, description) in _NONB_31.items():
            out += stanza([
                f"track nonb_{n}",
                "parent nonb on",
                f"shortLabel {label}",
                f"longLabel Non-B DNA: {description}",
                "type bigBed 6",
                f"bigDataUrl nonb_{n}.bb",
                "color 0,0,0",
            ], indent="    ")

        if sp == "Hsapien":
            out += stanza([
                "track g4_cuttag",
                "compositeTrack on",
                "shortLabel G4 CUT&Tag",
                "longLabel Published G4 CUT&Tag in HEK293T, mapped with the PDAL-Seq pipeline (CPM)",
                "type bigWig",
                "visibility hide",
                "autoScale on",
                "maxHeightPixels 100:40:8",
            ])
            for _, row in _G4_CUTTAG_31.iterrows():
                lab = str(row["PI"]).strip('"').split(",")[0]
                long_label = f"G4 CUT&Tag, HEK293T, {row['Run']} ({lab})"[:80]
                out += stanza([
                    f"track g4_cuttag_{row['Run']}",
                    "parent g4_cuttag on",
                    f"shortLabel {row['Run']}",
                    f"longLabel {long_label}",
                    "type bigWig",
                    f"bigDataUrl g4_cuttag_{row['Run']}.bw",
                    "color 0,0,0",
                ], indent="    ")

        with open(output[0], "w") as fh:
            fh.write(out.rstrip("\n") + "\n")
