import sys
import pandas as pd

def block_to_df(block_lines, i):
    """
    The block as a tidy formatted dataframe
    """
    sequence_lines = block_lines[1:-1]
    
    block = []
    species = []
    chrom = []
    start = []
    
    for s in sequence_lines:
        block.append(i)
        species.append(s.split("\t")[1].split(".")[0])
        chrom.append(s.split("\t")[1].split(".")[1])
        start.append(int(s.split("\t")[2]))
        
    data = {'block': block, 'species': species, 'chr': chrom, 'start': start}
    df = pd.DataFrame(data)
    
    return df

maf = sys.argv[1]
out = sys.argv[2]

block_open = False
block_lines = []
l_df = []

i = 1
with open(maf) as f:
    for line in f:
        LINE = line.strip()
        if LINE == "a":
            block_open = True
        if block_open:
            block_lines.append(LINE)
            if LINE == "":
                block_open = False
                l_df.append(block_to_df(block_lines, i))
                block_lines = [] #Delete the blocklines so we can start a new block
                i = i + 1

df = pd.concat(l_df)

df.to_csv(out, sep='\t')