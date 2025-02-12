import os
from Bio import AlignIO, SeqIO
import re

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
