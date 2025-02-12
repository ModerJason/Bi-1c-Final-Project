import os
from Bio import AlignIO, SeqIO
import re
import matplotlib.pyplot as plt

# File paths. The first is given in Stockholm form, while the second is given in FASTA form, which is desired.
sto_file = "pf00014/data/pf00014.sto"
fasta_file = "pf00014/data/pf00014.fasta"

# Check if the FASTA file already exists. If it already exists, then this line is ignored
if not os.path.exists(fasta_file):

    # Read the Stockholm alignment
    stockholm = AlignIO.read(sto_file, "stockholm")

    # Write to FASTA format
    with open(fasta_file, "w") as output_handle:
        AlignIO.write(stockholm, output_handle, "fasta")

fasta_file = "pf00014/data/pf00014.fasta"
fasta_filtered = "pf00014/data/filtered_pf00014.fasta"

# Check to see if filtered FASTA file already exists
if not os.path.exists(fasta_filtered):
    # Find BPT1_BOVIN sequence in pf00014.fasta, which is the guiding sequence for residues
    bovin_record = None
    for record in SeqIO.parse(fasta_file, "fasta"):
        if "BPT1_BOVIN" in record.id:  # Adjust if ID format differs
            bovin_record = record
            break

    # Identify non-gapped positions in BPT1_BOVIN
    non_gapped_positions = [i for i, residue in enumerate(bovin_record.seq) if residue != "-"]

    # Open output file for writing
    with open(fasta_filtered, "w") as out_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Extract only the non-gapped positions from each sequence
            filtered_sequence = "".join([record.seq[i] for i in non_gapped_positions])

            # Write the filtered sequence to the output file
            new_record = record[:]
            new_record.seq = filtered_sequence
            SeqIO.write(new_record, out_file, "fasta")

fasta_final = "pf00014/data/final_pf00014.fasta"

# Check to see if the filtered FASTA file already exists
if not os.path.exists(fasta_final):

    # Define the regex pattern for 7+ consecutive dashes (gaps)
    gap_pattern = re.compile(r"-{7,}")

    # Open output file for writing
    with open(fasta_final, "w") as out_file:
        # Iterate over each sequence in the input FASTA file
        for record in SeqIO.parse(fasta_filtered, "fasta"):
            # Convert sequence to uppercase and check for 7+ consecutive dashes
            sequence_upper = str(record.seq).upper()
            # Remove sequences with gaps or have an X or B (which are invalid in FASTA)
            if not gap_pattern.search(sequence_upper) and "X" not in sequence_upper and "B" not in sequence_upper:
                record.seq = sequence_upper  # Capitalize sequence so that these sequences are readable by adabmDCA
                SeqIO.write(record, out_file, "fasta")

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

# Plot the energies of all sequences in the PF00014 domain as well as the energy changed induced by point mutations
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

energies_fasta = "pf00014/bmDCA/energies.fasta"
sequence_ids = []
energies = []

# Read in the energies for each sequence in the PF00014 domain and plot frequencies in a histogram
with open(energies_fasta, "r") as file:
    for line in file:
        if line.startswith(">"):  # Header line
            parts = line.strip().split("|")  # Split by '|' to isolate energy value
            seq_id = parts[0][1:].split(" ")[0]  # Remove '>' and get sequence ID
            energy = float(parts[1].split(":")[1].strip())  # Extract energy value

            sequence_ids.append(seq_id)
            energies.append(energy)

axes[0].hist(energies, bins=100, color='blue', edgecolor='black', alpha=0.7)
axes[0].set_xlabel("Energy")
axes[0].set_ylabel("Frequency")
axes[0].set_title("bmDCA Energy Frequencies (PF00014)")
axes[0].grid(axis="y", linestyle="--", alpha=0.6)

mutations_fasta = "pf00014/bmDCA/mutation.fasta"
mutation_ids = []
delta_energies = []

# Read in energy changes caused by point-wise mutations and plot frequencies of changes in energy as a result
with open(mutations_fasta, "r") as file:
    for line in file:
        if line.startswith(">"):  # Header line
            parts = line.strip().split("|")  # Split by '|' to isolate energy value
            mutation_id = parts[0][1:].split(" ")[0]  # Remove '>' and get sequence ID
            delta_energy = float(parts[1].split(":")[1].strip())  # Extract energy value

            mutation_ids.append(mutation_id)
            delta_energies.append(delta_energy)

axes[1].hist(delta_energies, bins=100, color='blue', edgecolor='black', alpha=0.7)
axes[1].set_xlabel(r"$\Delta E$")
axes[1].set_ylabel("Frequency")
axes[1].set_title("Change in Energy Induced by Point Mutation of BPT1_BOVIN")
axes[1].grid(axis="y", linestyle="--", alpha=0.6)
plt.show()