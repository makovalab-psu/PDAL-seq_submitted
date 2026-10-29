import subprocess
import os
import csv

query_files    = [f'data/27_Updated_Human_HMM_functional_analysis/state_{i}.bed'
                  for i in range(0, 10)]
database_files = [f'data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_{i}.bed'
                  for i in range(0, 8)]

bash_script = 'workflow/scripts/09_Human_HMM/intersect_state_beds.sh'
os.chmod(bash_script, 0o755)

with open('data/27_Updated_Human_HMM_functional_analysis/intersect_results.csv',
          mode='w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(["query", "database", "total", "intersect"])

    for query in query_files:
        for database in database_files:
            try:
                result = subprocess.run(
                    [bash_script, query, database],
                    capture_output=True, text=True, check=True
                )
                total, intersect = result.stdout.strip().split(',')
                csv_writer.writerow([query, database, total, intersect])
            except subprocess.CalledProcessError as e:
                print(f'Error processing {query} and {database}: {e.stderr.strip()}')

print("Results written to data/27_Updated_Human_HMM_functional_analysis/intersect_results.csv")
