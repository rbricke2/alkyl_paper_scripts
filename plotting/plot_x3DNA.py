"""
    Author:          Rachel Bricker
    Year:            2025
"""

"""
   usage: python3 plot_x3DNA.py
      1. list of model names for legend (must be parallel w.r.t the arguments i to i+n)
      i. path to directories that contain arguments (2.) and (3.) that you want to plot

   example: python3 plot_x3DNA.py \
            "(a),(b),(c),(d),(e)" \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/dsDNA1/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/dsDNA3/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/dsDNA2/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/dsDNA4/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/dsDNA5/
            
            python3 plot_x3DNA.py \
            "(a)" \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/movies_gifs_pictures/dsDNA1_CHARMM/x3DNA/
"""

import sys
import numpy as np 
from functions_for_plots import *

# command line input
input_list = sys.argv[1].replace("\\n", "\n").replace("\\t", "\t").replace("\\(", "(").replace("\\)", ")")
legend     = input_list.split(',')
paths      = list(sys.argv[2:])

# add trailing forward slash to directory path if necessary
for path in range(len(paths)):
    path_split = paths[path].split("/")
    if path_split[-1] != "":
        paths[path] = paths[path] + "/"

################################################################################################
#
# FUNCTIONS
#
################################################################################################

def get_data(file_name):
    data = []
    time = []
    for i in range(1,3):
        print(i)
        n_skipped = 0
        skipping  = False # skip first time step of every file except the very first file (i.e. ignore repetitive data)
        with open(file_name+"_"+str(i)+".dat", "r") as f:
            for line in f:
                if len(line.strip()) != 0 : # make sure line isn't empty
                    line_split = line.split()
                    if line[0] in ['#']: # title, labels, and comments in file
                        print(line_split)
                        if line_split[1] == 'Time':
                            if n_skipped == 0 and i != 1:
                                skipping = True
                                n_skipped += 1
                            else:
                                time.append(int(float(line_split[3])))
                        else:
                            skipping = False
                        continue
                    if not skipping:
                        data.append([float(x) for x in line_split])
    return time, data

################################################################################################
#
# MAIN PROGRAM
#
################################################################################################

def main():
    time, data = get_data(paths[0]+"L-BPS")

if __name__ == "__main__": 
    main()
            