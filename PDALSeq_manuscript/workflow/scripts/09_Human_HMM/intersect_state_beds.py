import subprocess
import os
import csv

# Define the paths to your query and database files
query_files = [f'data/09_Human_HMM/04_Multivariate_guassian/10/state_{i}.bed' for i in range(0, 10)]  # Replace with your actual query file names
database_files = [f'data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_{i}.bed' for i in range(0, 8)]  # Replace with your actual database file names

# Path to your bash script
bash_script = 'workflow/scripts/09_Human_HMM/intersect_state_beds.sh'

# Ensure the bash script is executable
os.chmod(bash_script, 0o755)

# Open a CSV file for writing
with open('data/09_Human_HMM/04_Multivariate_guassian/10/intersect_results.csv', mode='w', newline='') as csvfile:
    # Create a CSV writer object
    csv_writer = csv.writer(csvfile)
    
    # Write the header row
    csv_writer.writerow(["query", "database", "total", "intersect"])

    # Iterate through each combination of query and database files
    for query in query_files:
        for database in database_files:
            # Call the bash script with the query and database files
            try:
                result = subprocess.run([bash_script, query, database], capture_output=True, text=True, check=True)
                # Parse the output
                total, intersect = result.stdout.strip().split(',')
                # Write the results to the CSV file
                csv_writer.writerow([query, database, total, intersect])
            except subprocess.CalledProcessError as e:
                print(f'Error processing {query} and {database}: {e.stderr.strip()}')

print("Results have been written to intersect_results.csv")

