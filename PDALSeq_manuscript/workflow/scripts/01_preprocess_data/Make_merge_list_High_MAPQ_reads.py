#Makes a list of bam files to merge with samtools merge
#sys.argv[1] The file index containing metadata for each sequencing sample
#sys.arv[2] The experiment that we want
#sys.arv[3] The output text file for the -b flag of samtools merge

import pandas as pd
import sys

print(sys.argv[1])
print(sys.argv[2])
print(sys.argv[3])

df = pd.read_csv(sys.argv[1], sep='\t')
filtered_df = df[df['Sequencing_sample'] != 'Ssyndactylus_Karenina_LCL_PDALseq_js3008h_20241218_SSY-Ka-20mM-2-3_S3_R1_001']
filtered_df = filtered_df[filtered_df['Experiment'] == sys.argv[2]]
Sequencing_samples = "PDALseq_data/Mapped_deduplicated_high_MapQ_bam/" + filtered_df["Sequencing_sample"] + ".bam"

with open(sys.argv[3], 'w+') as f:
    # write elements of list
    for items in Sequencing_samples:
        f.write('%s\n' %items)
    
    print("File written successfully")


