####################################################################
# Step 24 — HMM Browser Tracks
#
# Combines the 8 per-state BED files from the hominoidae 8-state
# Phylo-HMGP model into a single BED9 browser track with state names
# and itemRgb colors for IGV and the UCSC Genome Browser.
#
# Input BEDs are treated as ancient() so step 07 is never re-triggered.
# Output: results/24_HMM_Browser_tracks/hominoidae/hsapien.bed
####################################################################

_STATES_24 = ["0", "1", "2", "3", "4", "5", "6", "7"]

# Maps old state integer (string) → (browser track name, RGB color)
_STATE_INFO_24 = {
    "4": ("8-origin_4",   "244,166,55"),
    "7": ("7-recomb_7",   "219,88,41"),
    "3": ("6-gene_3",     "136,45,113"),
    "1": ("5-variable_1", "182,219,255"),
    "6": ("4-variable_6", "182,219,255"),
    "0": ("3-low_0",      "56,99,80"),
    "5": ("2-low_5",      "56,99,80"),
    "2": ("1-low_2",      "56,99,80"),
}


rule HMM_Browser_track_hominoidae:
    input:
        ancient(expand(
            "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_{state}.bed",
            state=_STATES_24,
        ))
    output:
        "results/24_HMM_Browser_tracks/hominoidae/hsapien.bed"
    log:
        "logs/24_HMM_Browser_tracks/hominoidae/hsapien.txt"
    run:
        import os, re

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)

        rows = []
        with open(log[0], "w") as logf:
            for bed_path in input:
                state = re.search(r"state_(\d+)\.bed$", bed_path).group(1)
                name, color = _STATE_INFO_24[state]
                n = 0
                with open(bed_path) as fh:
                    for line in fh:
                        if line.startswith("track") or line.startswith("#"):
                            continue
                        parts = line.rstrip("\n").split("\t")
                        chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                        rows.append((chrom, start, end, name, ".", ".", ".", ".", color))
                        n += 1
                logf.write(f"{bed_path}: {n} regions → name={name} color={color}\n")

        rows.sort(key=lambda r: (r[0], r[1]))

        with open(output[0], "w") as out:
            out.write(
                'track name="Ape_HMM_States" '
                'description="Ape HMM States (hominoidae, 8 states)" '
                'visibility=2 itemRgb="On"\n'
            )
            for r in rows:
                out.write("\t".join(str(x) for x in r) + "\n")
