import sys
import pandas as pd

def similar(seq1, seq2):
    return sum([1 for a, b in zip(seq1, seq2) if a == b])

def filter_block(block_lines):
    """
    Returns an alignment block with only one representative per species
    """
    sequence_lines = block_lines[1:-1]
    ref_line = sequence_lines[0]
    Sequence = []
    Species = []
    Score = []
    for s in sequence_lines:
        Species.append(s.split("\t")[1].split(".")[0])
        Sequence_s = s.split("\t")[6].upper()
        Sequence.append(Sequence_s)
        Score.append(similar(Sequence[0], Sequence_s))
    data = {'Species': Species, 'Sequence': Sequence, 'Score': Score, 'lines': sequence_lines}
    df = pd.DataFrame(data)
    df = df[df['Species'] != Species[0]]
    Species_out = list(set(df['Species']))
    out_lines = []
    out_lines.append("a")
    out_lines.append(ref_line)
    for s in Species_out:
        df_s = df[df['Species'] == s].sort_values("Score", ascending=False)
        out_lines.append(list(df_s['lines'])[0])
    out_lines.append("")
    return out_lines



maf = sys.argv[1]
out=sys.argv[2]


header = ["##maf version=1 scoring=N/A", ""]
with open(out, 'w+') as f:
    # write elements of list
    for items in header:
        f.write('%s\n' %items)

#Prune the alignment by block size and precence of all species

block_open = False
block_lines = []

with open(maf) as f:
    for line in f:
        LINE = line.strip()
        if LINE == "a":
            block_open = True
        if block_open:
            block_lines.append(LINE)
            if LINE == "":
                block_open = False
                new_block = filter_block(block_lines)
                with open(out, 'a') as f:
                        # write elements of list
                        for items in new_block:
                            f.write('%s\n' %items)
                block_lines = []

