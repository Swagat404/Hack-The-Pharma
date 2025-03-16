from sklearn.neighbors import KNeighborsRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import numpy as np
import statistics as stats
from sklearn.metrics import mean_absolute_error


X_merged_matrix = np.load("data_merged.npy", allow_pickle=True)
print(X_merged_matrix)  # (26520, 12), 26520 cleaned kinase-inhibitor pairings

# processed data files are in Google Colab or uploaded by other teammates

# Train model using a comprehensive set of features about the inhibitor and kinases:
# inhibitor fingerprints, kinase amino acid sequence reduced dimension embeddings,
# and additional kinase binary information (i.e. mutation status and kinase group membership)

y_cleaned = y_data_cleaned["Kd (nM)"].values

# print(X_comprehensive.shape) # same number of sample observations as what we want use to predict
# print(y_cleaned.shape)

# Should use same train-test data subsets as XG-Boost to ensure a fair comparison of the models
# and reproducibility
# Split availabile data
X_train, X_test, y_train, y_test = train_test_split(X_merged_matrix, y_cleaned,
                                                    test_size=0.2, random_state=42)

# K-Nearest Neighbors Regression
model = KNeighborsRegressor(n_neighbors = 5)
model.fit(X_train, y_train)
y_pred = model.predict(X_test)

## Evaluate the model
# MSE
mse_baseline = np.var(y_test)
rel_mse = mean_squared_error(y_test, y_pred) / mse_baseline

# MAE
  # Naive model assumes that the previous value is the prediction for the 
  # current time step
naive_forecast = np.abs(np.diff(y_test))
mae_naive = np.mean(naive_forecast)
rel_mae = mean_absolute_error(y_test, y_pred) / mae_naive
print(f"Relative MSE: {rel_mse}, Relative MAE: {rel_mae}")

# Calculate selectivity scores
n_kinases = 442
def selectivity_scores(matrix_y_pred, threshold = 3000):
  #### based on mapping, need to convert back to K_d array with
  #### kinases as rows and inhibitors as columns
  selectivity_scores = (matrix_y_pred < threshold).sum(axis=1) / n_kinases
  return(selectivity_scores)



# Calculate selectivity scores

# mapped predictions vector into matrix with kinase rows and inhibitor columns
matrix_pred_rd = np.load("/content/y_pred_matrix_baseline1.npy", allow_pickle=True)
# actually don't need the true R_d in matrix format because
# already given like that in Excel and already given S_inhibitor

print(matrix_pred_rd.shape) # (442, 12)

n_kinases = 442
def selectivity_scores(matrix_y_pred, threshold = 3000):
  matrix_y_pred = np.where(matrix_y_pred == None, np.nan, matrix_pred_rd).astype(float)
  matrix_y_pred = np.nan_to_num(matrix_pred_rd, nan=999999)  # Replace NaN with a high value so ignore it

  selectivity_scores = np.less(matrix_y_pred, threshold).sum(axis=1) / n_kinases
  return(selectivity_scores)

inhib_selectivity_3000 = selectivity_scores(mod_pred_rd)
inhib_selectivity_300 = selectivity_scores(mod_pred_rd, threshold = 300)

print(inhib_selectivity_3000)
print(inhib_selectivity_300)