import subprocess
import os
import csv
import sys
import pandas as pd

NONB = ["APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z"]
SPECIES = ["Hsapien", "Ptroglodytes","Ppaniscus","Ggorilla","Pabelii", "Ppygmaeus","Ssyndactylus"]
df_subclass = pd.read_csv("data/17_PDALseq_reads_in_CenSat/00_scratch/Filtered_subclasses_for_analysis.csv", sep=',') 
SUBCLASS = df_subclass['file_path'].tolist()

output_path = sys.argv[1]

# Path to your bash script
bash_script = 'workflow/scripts/17_PDALseq_reads_in_CenSat/intersect_censat_nonB.sh'

# Ensure the bash script is executable
os.chmod(bash_script, 0o755)

# Open a CSV file for writing
with open(output_path, mode='w', newline='') as csvfile:
    # Create a CSV writer object
    csv_writer = csv.writer(csvfile, delimiter='\t')
    # Write the header row
    csv_writer.writerow(["species", "censat_subclass", "nonb", "total", "intersect"])
    #loop through all combinations of nonb, species, and censats
    for nonb in NONB:
        for species in SPECIES:
            for censat in SUBCLASS:
                #Find censat file
                query = "data/17_PDALseq_reads_in_CenSat/03_censat_subclasses/"+species+"/"+censat+".bed"
                #Find nonb file
                database = "data/05_nonB_motif_cluster_analysis/01_Add_nonB_name/"+species+"/"+nonb+".bed"
                try:
                    result = subprocess.run([bash_script, query, database], capture_output=True, text=True, check=True)
                    # Parse the output
                    total, intersect = result.stdout.strip().split(',')
                    # Write the results to the CSV file
                    csv_writer.writerow([species, censat, nonb, total, intersect])
                    print(species+" "+censat+" "+nonb+" "+total+" "+intersect)
                except subprocess.CalledProcessError as e:
                    print(f'Error processing {query} and {database}: {e.stderr.strip()}')

print("Results have been written to intersect_results.csv")
