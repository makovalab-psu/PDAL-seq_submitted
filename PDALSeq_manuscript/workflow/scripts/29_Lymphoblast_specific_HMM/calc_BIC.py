import pandas as pd
import numpy as np
from hmmlearn import vhmm
import joblib
import sys

HMM_state = ["2", "4", "6", "8", "10", "12", "14", "16", "18", "20"]

model_path = sys.argv[1]
data_filename = sys.argv[2]
df = pd.read_csv(data_filename)

bic_values = []
log_l = []
params = []
log_obs = []

for x in HMM_state:
	print(x)
	model = joblib.load(model_path+x+"/GHMM.pkl")
	data_array = df.values
	log_likelihood = model.score(data_array)
	n_components = model.n_components
	n_features = model.n_features
	n_means = n_components * n_features
	n_covariances = n_components * (n_features * (n_features + 1)) // 2
	n_transitions = n_components * (n_components - 1)
	n_params = n_means + n_covariances + n_transitions
	n_observations = data_array.shape[0]
	bic = -2 * log_likelihood + n_params * np.log(n_observations)
	bic_values.append(bic)
	log_l.append(log_likelihood)
	params.append(n_params)
	log_obs.append(np.log(n_observations))

data_out = {"States":HMM_state, "BIC":bic_values, "log_l":log_l, "params":params, "log_obs":log_obs}
df_out = pd.DataFrame(data_out)
df_out.to_csv(sys.argv[3], sep=',', encoding='utf-8')
