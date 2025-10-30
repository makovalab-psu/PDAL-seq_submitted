import pandas as pd
import numpy as np
from hmmlearn import vhmm
import joblib
import sys

def main(model_file, transition_output):
    # Load model
    try:
        model = joblib.load(model_file)
    except Exception as e:
        print(f"Error loading model: {e}")
        return

    # Extract transition probability matrix
    transition_matrix = model.transmat_  # Shape: (n_components, n_components)
    
    # Create a tidy DataFrame
    n_components = model.n_components
    tidy_data = []

    for i in range(n_components):
        for j in range(n_components):
            tidy_data.append({
                'from_state': f'State_{i}',
                'to_state': f'State_{j}',
                'transition_probability': transition_matrix[i, j]
            })

    tidy_df = pd.DataFrame(tidy_data)

    # Prepare the transition data and write it to CSV
    tidy_df.to_csv(transition_output, sep=',', encoding='utf-8', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <model_file> <transition_output>")
    else:
        main(sys.argv[1], sys.argv[2])

