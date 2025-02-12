import matplotlib.pyplot as plt

threshold = 0.2
files = {
    "bmDCA": "pf00014/bmDCA/frobenius.txt",
    "eaDCA": "pf00014/eaDCA/eaDCA_frobenius.txt",
    "edDCA": "pf00014/edDCA/edDCA_frobenius.txt"
}

# Function to extract valid (i, j) pairs from a file
def extract_valid_pairs(filename, threshold):
    i_vals, j_vals = [], []
    with open(filename, "r") as file:
        for line in file:
            i, j, value = map(float, line.strip().split(","))
            if value > threshold and i < j:
                i_vals.append(int(i))
                j_vals.append(int(j))
    return i_vals, j_vals

# Extract data
bmDCA_i_vals, bmDCA_j_vals = extract_valid_pairs(files["bmDCA"], threshold)
eaDCA_i_vals, eaDCA_j_vals = extract_valid_pairs(files["eaDCA"], threshold)
edDCA_i_vals, edDCA_j_vals = extract_valid_pairs(files["edDCA"], threshold)

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Plot bmDCA
axes[0].scatter(bmDCA_i_vals, bmDCA_j_vals, color='purple', marker='d', s=10, label="bmDCA")
axes[0].set_title("bmDCA Predicted Contact Map")
axes[0].set_xlabel(r"$i$")
axes[0].set_ylabel(r"$j$")
axes[0].spines['bottom'].set_position('zero')

# Plot eaDCA
axes[1].scatter(eaDCA_i_vals, eaDCA_j_vals, color='blue', marker='s', s=5, label="eaDCA")
axes[1].set_title("eaDCA Predicted Contact Map")
axes[1].set_xlabel(r"$i$")
axes[1].set_ylabel(r"$j$")
axes[1].spines['left'].set_position('zero')
axes[1].spines['bottom'].set_position('zero')

# Plot edDCA
axes[2].scatter(edDCA_i_vals, edDCA_j_vals, color='green', s=5, label="edDCA")
axes[2].set_title("edDCA Predicted Contact Map")
axes[2].set_xlabel(r"$i$")
axes[2].set_ylabel(r"$j$")
axes[2].spines['left'].set_position('zero')
axes[2].spines['bottom'].set_position('zero')

plt.show()