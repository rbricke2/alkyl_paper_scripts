
# Overview 

Script used to construct the `.rtp` files used for CHARMM36 simulations that define the residue topology of the non-standard nucleotide fragments.

# Script

* `make_rtp_file.py`: Creates the `.rtp` files that define the residue topology of a central thymine, cytosine, guanine, or adenine fragment containing an ethyl- or decyl-phosphate group.

# Directories

* `rtp_mol2_files/`: Includes the `.mol2` files used as input for CGenFF and the `.rtp` files outputted by CGenFF (located in the GROMACS format file). The input `.mol2` files are needed to recover the atom names.
