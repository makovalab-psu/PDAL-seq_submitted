import sys
import pandas as pd

def similar(seq1, seq2):
    return sum([1 for a, b in zip(seq1, seq2) if a == b])/len(seq1)

def min_score(block_lines, taxon):
    """
    Returns the minimun fraction bit-score match to the reference for each sequence in the block
    """
    sequence_lines = block_lines[1:-1]
    Sequence = []
    Score = []
    Species = []
    for s in sequence_lines:
        Sequence_s = s.split("\t")[6].upper()
        Sequence.append(Sequence_s)
        Score.append(similar(Sequence[0], Sequence_s))
        Species.append(s.split("\t")[1].split(".")[0])
        # Find indices of elements in Species that are in taxon
    indices = [index for index, species in enumerate(Species) if species in taxon]
    # Extract elements at those indices
    result_Score = [Score[i] for i in indices]
    return min(Score)

maf = sys.argv[1]
out = sys.argv[2]
MIN = float(sys.argv[3])
TAXON = sys.argv[4]

#hominoidae
if (TAXON == "hominoidae"):
    taxon = ["hs1", "GCA_029281585", "GCA_028858775", "GCA_029289425", "GCA_028885655", "GCA_028885625", "GCA_028878055"]
#hominidae
if (TAXON == "hominidae"):
    taxon = ["hs1", "GCA_029281585", "GCA_028858775", "GCA_029289425", "GCA_028885655", "GCA_028885625"]
#homininae
if (TAXON == "homininae"):
    taxon =  ["hs1", "GCA_029281585", "GCA_028858775", "GCA_029289425"]
#hominini
if (TAXON == "hominini"):
    taxon =  ["hs1", "GCA_028858775", "GCA_029289425"]

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
                write_block = True
                MIN_score = min_score(block_lines, taxon)
                if MIN_score <= MIN:
                    write_block = False
                if write_block: #Write the block if it passes filters
                    with open(out, 'a') as f:
                        # write elements of list
                        for items in block_lines:
                            f.write('%s\n' %items)
                    write_block = False
                block_lines = [] #Delete the blocklines so we can start a new block


