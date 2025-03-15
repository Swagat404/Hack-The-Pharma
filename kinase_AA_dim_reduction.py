import numpy as np


kinase_embedd = np.load("global_embeddings_max_length.npy")
# print(kinase_embedd.shape) # Original shape ~ (442, 15599)

# Compute full svd for comparison
U, S, Vt = np.linalg.svd(kinase_embedd)

cumulative_variance = np.cumsum(S**2) / np.sum(S**2) # cumulative up to each k

# Choose k for 95% variance retention
optimal_k = np.argmax(cumulative_variance >= 0.95) + 1
print(f"Optimal k for 95% variance retention: {optimal_k}")

# So let's choose projecting onto a 38-dimensional space

# If wanted to check explained variance for a specific k
def explained_variance(k):
    explained_variance = cumulative_variance[k-1]
    print(f"Explained Variance for k={k}: {explained_variance * 100:.2f}%")

svd = TruncatedSVD(n_components = optimal_k)  # Reduce to optimal k dimensions
embedd_reduced = svd.fit_transform(kinase_embedd)

print("Original shape:", kinase_embedd.shape)  # (442, 15599)
print("Reduced shape:", embedd_reduced.shape)  # (442, 38)


