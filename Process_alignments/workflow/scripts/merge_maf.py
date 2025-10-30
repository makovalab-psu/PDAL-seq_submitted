import sys
import pandas as pd

maf_dir = sys.argv[1]
out = sys.argv[2]

CHR = pd.read_csv("resources/chr_sizes/Hsapien.chrom.sizes", sep='\t', header=None)[0].tolist()

header = ["##maf version=1 scoring=N/A", ""]

with open(out, 'w+') as f:
        # write elements of list
        for items in header:
            f.write('%s\n' %items)


for chr in CHR:
    maf = maf_dir+"/"+chr+".maf"
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
                with open(out, 'a') as f:
                        # write elements of list
                        for items in block_lines:
                            f.write('%s\n' %items)
                block_lines = [] #Delete the blocklines so we can start a new block


