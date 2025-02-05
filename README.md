# Bi-1c-Final-Project
Explorations in adabmDCA 2.0, based on the paper https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04441-9.

In this project, we explore the new adabmDCA 2.0 package (https://github.com/spqb/adabmDCA), which is still being actively updated. We focus on the PF00014 (Kunitz) domain. We have done multiple things in this project.

- Test and plot protein contact maps for three different Boltzmann machine Direct Coupling Analysis (DCA) algorithms:
  - bmDCA (Boltzmann machine), a fully-connected DCA model
  We also test two sparser Boltzmann machine DCA methods, which is well-motivated due to the relatively sparse interaction network of amino acid sequences.
  - eaDCA (element activation), which starts with an empty set of connections, and systematically adds couplings every epoch
  - edDCA (element decimation), which starts from a fully-connected graph and removes couplings which produce the smallest peturbations on the probability distributions every epoch
- Calculate and plot energies of each sequence in the PF00014 domain as determined by the trained bmDCA model, and determine the effect on energy caused by point mutations in the BPT1_BOVIN 39-91 sequence. The BPT1_BOVIN 39-91 sequence is chosen because this protein specificies the relevant residues in the FASTA file.
