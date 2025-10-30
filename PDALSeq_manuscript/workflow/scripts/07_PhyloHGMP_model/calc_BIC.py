import pandas as pd
import numpy as np
from hmmlearn import vhmm
import joblib
import sys

#Make a list of state numbers
HMM_state = ["2", "4", "6", "8", "10", "12", "14", "16", "18", "20", "25", "30"]

# Load the model from the specified file
model_path = sys.argv[1]  # Model filename from command line argument


# Load the data for which to calculate BIC
data_filename = sys.argv[2]  # Data filename from command line argument
df = pd.read_csv(data_filename)

# Initialize a list to store tuples of (model_filename, BIC)
bic_values = []
log_l = []
params = []
log_obs = []

for x in HMM_state:
	print(x)
	model = joblib.load(model_path+x+"/GHMM.pkl")
	# Convert the DataFrame to a numpy array
	data_array = df.values
	# Calculate the log likelihood of the data given the model
	log_likelihood = model.score(data_array)
	# Number of parameters in the model for multivariate Gaussian HMM
	n_components = model.n_components
	n_features = model.n_features
	# Mean parameters
	n_means = n_components * n_features
	# Covariance parameters (unique elements in the covariance matrix)
	n_covariances = n_components * (n_features * (n_features + 1)) // 2
	# Transition parameters
	n_transitions = n_components * (n_components - 1)
	# Total number of parameters
	n_params = n_means + n_covariances + n_transitions
	# Number of observations
	n_observations = data_array.shape[0]
	# Calculate BIC
	bic = -2 * log_likelihood + n_params * np.log(n_observations)
	# Save the BIC value
	bic_values.append(bic)
	log_l.append(log_likelihood)
	params.append(n_params)
	log_obs.append(np.log(n_observations))

#Prepare the BIC data and write it
data_out = {"States":HMM_state, "BIC":bic_values, "log_l":log_l, "params":params, "log_obs":log_obs}
df_out = pd.DataFrame(data_out)

df_out.to_csv(sys.argv[3], sep=',', encoding='utf-8')

