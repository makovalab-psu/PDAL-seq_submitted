import sys

def sort_a_by_indecies_in_b(a, b):
    # Create a list of indices based on the sorted order of b
    sorted_indices = sorted(range(len(b)), key=lambda x: b[x])
    # Use the sorted indices to reorder list a
    sorted_a = [a[i] for i in sorted_indices]
    return sorted_a

#test sort_a_by_indecies_in_b
a = ['apple', 'banana', 'cherry', 'date', "fish"]
b = [1, 2, 4, 3, 3]
#should return ['apple', 'banana', 'date', "fish", 'cherry']
sort_a_by_indecies_in_b(a, b)

maf = sys.argv[1]
out = sys.argv[2]

block_open = False
block_lines = []
l_block_lines = []
block_start = []

with open(maf) as f:
    for line in f:
        LINE = line.strip()
        if LINE == "a":
            block_open = True
        if block_open:
            block_lines.append(LINE)
            if LINE == "":
                block_open = False
                block_start.append(int(block_lines[1].split("\t")[2]))
                l_block_lines.append(block_lines)
                block_lines = [] #Delete the blocklines so we can start a new block

sorted_l_block_lines = sort_a_by_indecies_in_b(l_block_lines, block_start)


header = ["##maf version=1 scoring=N/A", ""]
with open(out, 'w+') as f:
    # write elements of list
    for items in header:
        f.write('%s\n' %items)    

for block_line in l_block_lines:
    with open(out, 'a') as f:
        # write elements of list
        for items in block_line:
            f.write('%s\n' %items)