import sys

####This bit was provided by chat gpt####

#I need a function to test if all species are present in each alignment block. This is a problem where we need to test if each string in a list is present as a substring in another list.
def all_strings_are_not_substrings(list1, list2):
    """
    Returns True if not every string in list1 is a substring of at least
    one string in list2. Otherwise, returns False.
    """
    for s in list1:
        if not any(s in t for t in list2):
            return True
    return False

#Test
list1 = ["cat", "dog"]
list2 = ["my cat is cute", "a big dog", "bird"]
all_strings_are_not_substrings(list1, list2)
#Further testing this chat gpt function
list3 = ["my bird is cute", "a big dog", "bird"]
all_strings_are_not_substrings(list1, list3)
list4 = ["mycatiscute", "abigdog", "bird"]
all_strings_are_not_substrings(list1, list4)
####Done with chat gpt####

#Get the system arguments
maf = sys.argv[1]
out=sys.argv[2]
Min_bp = int(sys.argv[3])

#Write the header
header = ["##maf version=1 scoring=N/A", ""]
with open(out, 'w+') as f:
    # write elements of list
    for items in header:
        f.write('%s\n' %items)

#Prune the alignment by block size and precence of all species
species = ["hs1", "GCA_028858775", "GCA_029289425"]
block_open = False
block_lines = []

with open(maf) as f:
    for line in f:
        LINE = line.strip() #Read one line from the file
        if LINE == "a":
            block_open = True #Start a alignment block
        if block_open:
            block_lines.append(LINE) #Start recording the the block
            if LINE == "": #End the block and write the block to a new file if it passes filters
                block_open = False
                write_block = True
                block_length = len(block_lines[1].split("\t")[6]) # determine the block length
                if block_length <=  Min_bp: #Determine if the block is long enough
                    write_block = False
                if all_strings_are_not_substrings(species, block_lines): #Determine if all species are represented in a block
                    write_block = False
                if write_block: #Write the block if it passes filters
                    with open(out, 'a') as f:
                        # write elements of list
                        for items in block_lines:
                            f.write('%s\n' %items)
                    write_block = False
                block_lines = [] #Delete the blocklines so we can start a new block
