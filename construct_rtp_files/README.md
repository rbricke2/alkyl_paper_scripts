
# Overview 

Script used to construct the `.rtp` files used for CHARMM36 GROMACS simulations that define the residue topology of the non-standard nucleotide fragments.

# Script

* `make_rtp_file.py`: Creates the `.rtp` files that define the residue topology of a central thymine, cytosine, guanine, or adenine fragment containing an ethyl- or decyl-phosphate group.

# Directories

* `rtp_mol2_files/`: Includes the `.mol2` files used as input for CGenFF and the `.rtp` files outputted by the <cite>[CGenFF server][1]</cite> (located in the GROMACS format file). The input `.mol2` files are needed to recover the atom names.

# Acknowlegdements

Thanks to Dr. Aleksei Aksimentiev and Dr. Christopher Maffeo for helping me with parameterizing a central fragment of a nucleotide containing an alkyl-phosphate group using CGenFF server.

[1]: https://doi.org/10.1002/jcc.21367
