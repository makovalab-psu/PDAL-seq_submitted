"""Write the ape HMM states for one species as BED9 with Fig. 3 labels and colors.

Usage: hmm_states_to_bed9.py GHMM_STATES_CSV SPECIES_ALIGNMENT_BED > out.bed

Row i of GHMM_states.csv is alignment block i, and every species' hominoidae
alignment BED lists the blocks in the same order (see Normalize_alignment_coverage.R
and compile_alignment_coverage.R in step 07), so the state of block i applies to
row i of each species' BED.
"""
import csv
import sys

# From Revised_figures/Figure_3_Ape_PDAL-Seq_MVGHMM.R:
#   State.levels = c("2", "5", "0", "6", "1", "3", "7", "4")
#   State.labels = rev(c("8-origin", "7-recomb", "6-gene", "5-variable", "4-variable", "3-low", "2-low", "1-low"))
#   scale_fill_manual(values = c("#386350", "#386350", "#386350", "#B6DBFF", "#B6DBFF", "#882D71", "#DB5829", "#F4A637"))
STATE_LEVELS = ["2", "5", "0", "6", "1", "3", "7", "4"]
STATE_LABELS = ["1-low", "2-low", "3-low", "4-variable", "5-variable", "6-gene", "7-recomb", "8-origin"]
STATE_COLORS = ["#386350", "#386350", "#386350", "#B6DBFF", "#B6DBFF", "#882D71", "#DB5829", "#F4A637"]


def hex_to_rgb(color):
    return ",".join(str(int(color[i:i + 2], 16)) for i in (1, 3, 5))


STATES = {s: (label, hex_to_rgb(color)) for s, label, color in zip(STATE_LEVELS, STATE_LABELS, STATE_COLORS)}

with open(sys.argv[1]) as fh:
    states = [str(int(float(row["hidden_state"]))) for row in csv.DictReader(fh)]

with open(sys.argv[2]) as fh:
    blocks = [line.rstrip("\n").split("\t")[:3] for line in fh if line.strip()]

if len(states) != len(blocks):
    sys.exit(f"{len(states)} states but {len(blocks)} alignment blocks in {sys.argv[2]}")

for (chrom, start, end), state in zip(blocks, states):
    label, rgb = STATES[state]
    sys.stdout.write(f"{chrom}\t{start}\t{end}\t{label}\t0\t.\t{start}\t{end}\t{rgb}\n")
