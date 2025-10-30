import sys
import pandas as pd

maf = sys.argv[1]
chr_names = sys.argv[2]
out = sys.argv[3]

CHR = pd.read_csv("resources/chr_sizes/Hsapien.chrom.sizes", sep='\t', header=None)[0].tolist()

header = ["##maf version=1 scoring=N/A", ""]

for chr in CHR:
    out_file = out+"/"+chr+".maf"
    with open(out_file, 'w+') as f:
        # write elements of list
        for items in header:
            f.write('%s\n' %items)

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
                chr = block_lines[1].split("\t")[1].split(".")[1]
                out_file = out+"/"+chr+".maf"
                with open(out_file, 'a') as f:
                        # write elements of list
                        for items in block_lines:
                            f.write('%s\n' %items)
                block_lines = [] #Delete the blocklines so we can start a new block


