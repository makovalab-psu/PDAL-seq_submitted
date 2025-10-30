import pandas as pd
import numpy as np
from hmmlearn import vhmm
import joblib
import sys


#Load model
model = joblib.load(sys.argv[1])

#extract n_features
# Number of parameters in the model for multivariate Gaussian HMM
n_components = model.n_components
n_features = model.n_features
x = n_components
	
# Extract mean parameters
means = model.means_  # Shape: (n_components, n_features)
# Convert means to a DataFrame
means_df = pd.DataFrame(means, columns=[f'Feature_{i+1}' for i in range(n_features)])
means_df['state'] = [f'{i}' for i in range(n_components)]
means_df['n_states'] = [f'{x}' for i in range(n_components)]

# Extract variance parameters
variances = np.array([np.diag(cov) for cov in model.covars_])  # Shape: (n_components, n_features)
variances_df = pd.DataFrame(variances, columns=[f'Variance_{i+1}' for i in range(n_features)])
variances_df['state'] = [f'{i}' for i in range(n_components)]
variances_df['n_states'] = [f'{x}' for i in range(n_components)]  # Assign the current state number

#Prepare the means data and write it
means_df.to_csv(sys.argv[3], sep=',', encoding='utf-8')

#Prepare the varience data and write it
variances_df.to_csv(sys.argv[4], sep=',', encoding='utf-8')

