import matplotlib.pyplot as plt

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