#!/usr/bin/env python3
"""
profile_from_matrix.py
Streaming mean profile from deeptools computeMatrix .mat.gz files.

Three modes:

  Accumulate (no flag):
    python profile_from_matrix.py <matrix.mat.gz> <state.json>
    Reads the matrix line-by-line and merges values into a running-sum state
    file.  Creates state.json on first call; updates it on subsequent calls.
    Peak RAM = O(N_bins × N_samples) — only the running sums, never the full
    feature matrix.

  Finalize (--finalize):
    python profile_from_matrix.py --finalize <state.json> <output.tsv>
    Divides accumulated sums by counts and writes the mean-profile TSV with
    columns: Bin Labels, forward, reverse.
    Used for BED3 (unstranded) features.

  Combine strands (--combine-strands):
    python profile_from_matrix.py --combine-strands <state_plus.json> <state_minus.json> <output.tsv>
    Merges two per-strand accumulation states into a single template/nontemplate
    profile TSV with columns: Bin Labels, template, nontemplate.
    Used for BED6 (strand-aware) features after separate + and − accumulation.

    Recoding rationale (G4/GQ example):
      + strand feature: G4 forms on + strand; displaced ssDNA is on − strand
                        → template = reverse bigwig (sums[1])
                        → nontemplate = forward bigwig (sums[0])
      − strand feature: deeptools has already flipped window 5'→3', so bin
                        positions are comparable to + strand features.
                        G4 forms on − strand; displaced ssDNA is on + strand
                        → template = forward bigwig (sums[0])
                        → nontemplate = reverse bigwig (sums[1])
    Sums and counts from both states are combined before dividing, so the
    resulting mean is weighted by feature count, not a mean of means.

Called from compute_strand_profile.sh.
"""

import gzip
import json
import math
import os
import sys


def _scalar(v, default):
    """Unwrap a deeptools list-valued header field (e.g. upstream=[1000]) to int."""
    if v is None:
        return default
    return v[0] if isinstance(v, list) else v


def accumulate(matrix_file, state_file):
    """Stream one .mat.gz and merge its values into the running state."""
    with gzip.open(matrix_file, "rt") as fh:
        # ── Parse JSON header (first line, prefixed '@') ──────────────────
        header_line = fh.readline()
        if not header_line.startswith("@"):
            sys.exit(f"ERROR: expected '@' header, got: {header_line[:80]!r}")
        params = json.loads(header_line[1:])

        sample_labels     = params["sample_labels"]
        upstream          = int(_scalar(params.get("upstream"),   1000))
        downstream        = int(_scalar(params.get("downstream"), 1000))
        bin_size          = int(_scalar(params.get("bin size"),   10))
        n_samples         = len(sample_labels)
        n_bins            = (upstream + downstream) // bin_size
        sample_boundaries = params.get(
            "sample_boundaries",
            list(range(0, n_samples * n_bins + 1, n_bins)),
        )

        # ── Load existing state or initialize ─────────────────────────────
        if os.path.exists(state_file):
            with open(state_file) as sf:
                state = json.load(sf)
            sums = state["sums"]
            cnts = state["cnts"]
        else:
            sums = [[0.0] * n_bins for _ in range(n_samples)]
            cnts = [[0]   * n_bins for _ in range(n_samples)]

        # ── Stream data rows (6 BED metadata columns, then coverage) ──────
        N_META = 6
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            data = cols[N_META:]
            for s in range(n_samples):
                start_col = sample_boundaries[s]
                end_col   = sample_boundaries[s + 1]
                for b_local, col_idx in enumerate(range(start_col, end_col)):
                    if col_idx >= len(data):
                        break
                    raw = data[col_idx]
                    if not raw:
                        continue
                    try:
                        v = float(raw)
                    except ValueError:
                        continue
                    if not math.isnan(v):
                        sums[s][b_local] += v
                        cnts[s][b_local] += 1

    # ── Persist updated state ──────────────────────────────────────────────
    with open(state_file, "w") as sf:
        json.dump(
            {
                "sample_labels": sample_labels,
                "upstream":      upstream,
                "downstream":    downstream,
                "bin_size":      bin_size,
                "n_bins":        n_bins,
                "n_samples":     n_samples,
                "sums":          sums,
                "cnts":          cnts,
            },
            sf,
        )


def finalize(state_file, output_tsv):
    """Divide accumulated sums by counts and write the mean-profile TSV."""
    with open(state_file) as sf:
        state = json.load(sf)

    sample_labels = state["sample_labels"]
    upstream      = state["upstream"]
    bin_size      = state["bin_size"]
    n_bins        = state["n_bins"]
    n_samples     = state["n_samples"]
    sums          = state["sums"]
    cnts          = state["cnts"]

    bin_centers = [
        int(-upstream + bin_size * b + bin_size // 2)
        for b in range(n_bins)
    ]

    with open(output_tsv, "w") as out:
        out.write("Bin Labels\t" + "\t".join(sample_labels) + "\n")
        for b in range(n_bins):
            row = [str(bin_centers[b])]
            for s in range(n_samples):
                row.append(
                    f"{sums[s][b] / cnts[s][b]:.6f}" if cnts[s][b] > 0 else "nan"
                )
            out.write("\t".join(row) + "\n")


def combine_strands(state_plus_file, state_minus_file, output_tsv):
    """
    Merge per-strand accumulation states into a single template/nontemplate TSV.

    sums[0] = forward bigwig, sums[1] = reverse bigwig (order matches -S flag in
    compute_strand_profile.sh: forward_bw then reverse_bw).

    + strand features: template = sums[1] (reverse), nontemplate = sums[0] (forward)
    − strand features: template = sums[0] (forward), nontemplate = sums[1] (reverse)
    """
    plus_exists  = os.path.exists(state_plus_file)
    minus_exists = os.path.exists(state_minus_file)

    if not plus_exists and not minus_exists:
        open(output_tsv, "w").close()
        return

    ref_file = state_plus_file if plus_exists else state_minus_file
    with open(ref_file) as sf:
        ref = json.load(sf)
    n_bins   = ref["n_bins"]
    upstream = ref["upstream"]
    bin_size = ref["bin_size"]

    template_sums    = [0.0] * n_bins
    template_cnts    = [0]   * n_bins
    nontemplate_sums = [0.0] * n_bins
    nontemplate_cnts = [0]   * n_bins

    if plus_exists:
        with open(state_plus_file) as sf:
            plus = json.load(sf)
        for b in range(n_bins):
            template_sums[b]    += plus["sums"][1][b]
            template_cnts[b]    += plus["cnts"][1][b]
            nontemplate_sums[b] += plus["sums"][0][b]
            nontemplate_cnts[b] += plus["cnts"][0][b]

    if minus_exists:
        with open(state_minus_file) as sf:
            minus = json.load(sf)
        for b in range(n_bins):
            template_sums[b]    += minus["sums"][0][b]
            template_cnts[b]    += minus["cnts"][0][b]
            nontemplate_sums[b] += minus["sums"][1][b]
            nontemplate_cnts[b] += minus["cnts"][1][b]

    bin_centers = [
        int(-upstream + bin_size * b + bin_size // 2)
        for b in range(n_bins)
    ]

    with open(output_tsv, "w") as out:
        out.write("Bin Labels\ttemplate\tnontemplate\n")
        for b in range(n_bins):
            t  = f"{template_sums[b]    / template_cnts[b]:.6f}"    if template_cnts[b]    > 0 else "nan"
            nt = f"{nontemplate_sums[b] / nontemplate_cnts[b]:.6f}" if nontemplate_cnts[b] > 0 else "nan"
            out.write(f"{bin_centers[b]}\t{t}\t{nt}\n")


def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print(__doc__)
        sys.exit(1)

    if args[0] == "--finalize":
        if len(args) < 3:
            sys.exit("Usage: profile_from_matrix.py --finalize <state.json> <output.tsv>")
        finalize(args[1], args[2])
    elif args[0] == "--combine-strands":
        if len(args) < 4:
            sys.exit("Usage: profile_from_matrix.py --combine-strands <state_plus.json> <state_minus.json> <output.tsv>")
        combine_strands(args[1], args[2], args[3])
    else:
        accumulate(args[0], args[1])


if __name__ == "__main__":
    main()
