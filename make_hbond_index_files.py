"""
    Author: Rachel Bricker
    Year:   2025
"""

"""
    Creates index files which will be used as input for GROMACS utilities 'angle' and 'distance'.
    
    THIS CODE REQUIRES:
        * ENERGY MINIMIZATION TO BE PERFORMED. HENCE, IT USES THE .gro FILE OUTPUTTED BY THE
          MDRUN THAT CARRIED OUT ENERGY MINIMIZATION
        * DNA IS DOUBLE STRANDED
"""

"""
   usage: python3 make_hbond_index_files.py
      1. [path to .gro file outputted after minimization]
      2. [number of base pairs in duplex]
      3. [output directory]
   
   example: python3 make_hbond_index_files.py /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/dsDNA1/em.gro 42 /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/dsDNA1/
"""

import sys
import numpy as np

# command line input
input_gro_file_path = str(sys.argv[1])
max_residue_id      = int(sys.argv[2])
output_dir          = str(sys.argv[3])

output_dir_split = output_dir.split("/")
if output_dir_split[-1] != "":
    output_dir = output_dir + "/"

################################################################################################
#
# FUNCTIONS
#
################################################################################################

def make_hbond_index_files(input_gro_file_path, max_residue_id, output_dir):
    # initialize dictionary which will track the atom id of the nitrogen atom in the nucleobase ring
    # in each residue that we will use to characterize hydrogen bonding
    atoms = {}

    ################################################
    # read input .gro file

    with open(input_gro_file_path, "r") as f:
        # initialize residue id tracker
        residue_id = 1
        
        # initialize residue name tracker
        residue_name = None
        
        # initialize ring type tracker
        purine = False

        for line in f:
            line_split = line.split()
            if len(line_split) > 1:    # make sure line isn't empty
                
                # check if we are at the next residue
                if (line_split[0] != residue_name) and (residue_name != None):
                    residue_id += 1        # increment residue id tracker
                    purine = False         # reset ring type tracker
                    if residue_id >= (max_residue_id+1):
                        break    # if we have iterated through the entire DNA duplex, then stop reading lines
                
                
                if line_split[0][0:len(str(residue_id))] == str(residue_id):
                    residue_name = line_split[0]    # get residue name
                    if line_split[1] == 'N9':
                        purine = True    # if N9 atom is in residue, then nucleobase is purine
                    elif (line_split[1] == 'N3') and (not purine):
                        atoms[residue_name] = line_split[2]    # N3 atom is the nitrogen atom in the nucleobase ring for pyrimidines, so record the atom ID
                    elif (line_split[1] == 'N1') and purine: 
                        atoms[residue_name] = line_split[2]    # N1 atom is the nitrogen atom in the nucleobase ring for purines, so record the atom ID

    ################################################
    # create index files
    
    # list of all recorded atom IDs
    all_ids = list(atoms.values())
    
    # list of all recorded residue names
    all_resi_names = list(atoms.keys())
    
    # compute midpoint (essentially the number of nucleotides/residues in each strand)
    mid_point = len(all_ids)//2
    
    # list of all recorded atom IDs in the first strand
    first_strand_ids = all_ids[:mid_point]
    
    # list of all recorded residue names in the first strand
    first_strand_resi_names = all_resi_names[:mid_point]
    
    # list of all recorded atom IDs in the second strand
    second_strand_ids = all_ids[mid_point:]
    
    # reverse list so that it is the reverse-complement counterpart of first_strand_ids
    second_strand_ids.reverse()
    
    
    file_name = "hbond_dist.ndx"    # file used as input for GROMACS utility 'distance'
    with open(output_dir + file_name, "w+") as output:
        output.write("[ H_BOND_FOR_DIST ]\n")
        for i in range(mid_point):
            # write the atom ids of the nitrogen atoms in each base pair
            output.write(first_strand_ids[i] + " " + second_strand_ids[i] + "\n")

    file_name = "hbond_angle.ndx"    # file used as input for GROMACS utility 'angle'
    with open(output_dir + file_name, "w+") as output:
        output.write("[ H_BOND_FOR_ANGLE ]\n")
        for i in range(mid_point):
            if 'G' in first_strand_resi_names[i] or 'T' in first_strand_resi_names[i]:
                # if the nucleobase of the first strand is guanine or thymine, then the first strand contains the donor/hydrogen;
                # thus, use the atom id of the nitrogen atom in the first strand to obtain the atom id of the hydrogen atom
                output.write(str(int(first_strand_ids[i])+1) + " " + first_strand_ids[i] + " " + second_strand_ids[i] + "\n")
            else:
                # else, the nucleobase of the second strand is guanine or thymine and the second strand contains the donor/hydrogen;
                # thus, use the atom id of the nitrogen atom in the second strand to obtain the atom id of the hydrogen atom
                output.write(str(int(second_strand_ids[i])+1) + " " + second_strand_ids[i] + " " + first_strand_ids[i] + "\n")

################################################################################################
#
# MAIN PROGRAM
#
################################################################################################

def main():
    make_hbond_index_files(input_gro_file_path, max_residue_id, output_dir)

if __name__ == "__main__": 
    main()
