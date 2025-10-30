import pandas as pd
import numpy as np
from hmmlearn import hmm
from hmmlearn import vhmm
import joblib
import sys

df = pd.read_csv(sys.argv[1])

# Convert the DataFrame to a numpy array
data_array = df.values

# Initialize the Gaussian HMM
n_states = int(sys.argv[2])  # Number of hidden states

model = vhmm.VariationalGaussianHMM(n_components = n_states, covariance_type='full', n_iter=10, verbose=True)

# Fit the model to the data
model.fit(data_array)

# Predict the hidden states for each observation
hidden_states = model.predict(data_array)

# Add the hidden states to the original DataFrame
df['hidden_state'] = hidden_states

# Display the first few rows of the DataFrame with hidden states
print(df.head())

# Save the DataFrame with hidden states to a CSV file
df.to_csv(sys.argv[3], index=False)

# Save the model to disk
joblib.dump(model, sys.argv[4])

