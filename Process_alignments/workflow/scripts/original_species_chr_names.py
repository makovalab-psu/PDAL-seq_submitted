import sys

maf = sys.argv[1]
species_out=sys.argv[2]
chr_out = sys.argv[3]

species = []
chr = []
with open(maf) as f:
    for line in f:
        LINE = line.strip()
        if LINE != "a" and LINE != "" and LINE != "##maf version=1 scoring=N/A":
            species_chr = LINE.split("\t")[1].split(".")
            if len(species_chr) == 2:
                species.append(species_chr[0])
                chr.append(species_chr[1])
            if len(species_chr) == 4:
                species.append(species_chr[0]+"."+species_chr[1])
                chr.append(species_chr[2]+"."+species_chr[3])


species_set = list(set(species))
chr_set = list(set(chr))

with open(species_out, 'w+') as f:
    # write elements of list
    for items in species_set:
        f.write('%s\n' %items)

with open(chr_out, 'w+') as f:
    # write elements of list
    for items in chr_set:
        f.write('%s\n' %items)