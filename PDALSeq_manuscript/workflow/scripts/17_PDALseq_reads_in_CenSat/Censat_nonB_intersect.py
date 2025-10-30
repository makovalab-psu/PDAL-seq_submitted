import subprocess
import os
import csv
import sys

NONB = ["APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z"]
SPECIES = ["Hsapien", "Ptroglodytes","Ppaniscus","Ggorilla","Pabelii", "Ppygmaeus","Ssyndactylus"]
CENSAT = ["active_hor", "bSat", "cenSat", "ct", "dhor", "GAP", "GAP,HSat3", "gSat", "hor", "HSat1A", "HSat1B", "HSat2", "HSat2_3", "HSat3", "mixedSuperFamily", "mon", "mon-hor", "rDNA", "subTerm", "subTerm_mon", "subTerm_mon-hor"]

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
    csv_writer.writerow(["species", "censat", "nonb", "total", "intersect"])
    #loop through all combinations of nonb, species, and censats
    for nonb in NONB:
        for species in SPECIES:
            for censat in CENSAT:
                #Find censat file
                query = "data/17_PDALseq_reads_in_CenSat/01_censat_classes/"+species+"/"+censat+".bed"
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

