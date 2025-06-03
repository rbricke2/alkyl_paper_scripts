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
import math
import copy
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

# this function was grabbed from stackoverflow:
# https://stackoverflow.com/questions/1628386/normalise-orientation-between-0-and-360
def normalize_range(old_angle : np.float64, start : float, end : float) -> np.float64:
    """
        Normalizes any number to an arbitrary range by assuming the range wraps around when going below the
        minimum or above the maximum

        Parameters:
            old_angle (numpy.float64) : the angle in the range [-180, 180]
            start     (float)         : the number that defines the first value in the new range
            end       (float)         : the number that defines the last value in the new range 

        Returns:
            new_angle (numpy.float64) : the angle in the range [`start`, `end`]
    """

    width         = end - start
    off_set_value = old_angle - start    # the old angle relative to zero
    new_angle     = (off_set_value - ( math.floor( off_set_value / width ) * width )) + start    # + start to reset back to start of original range

    return new_angle

def get_data(file_name):
    #Shift     Slide      Rise      Tilt      Roll     Twist
    dict_temp   = {"shift": [], "slide": [], "rise": [], "tilt": [], "roll": [], "twist": []}
    parameters  = list(dict_temp.keys())
    data        = {}
    time_step   = None
    for i in range(1,3):
        n_skipped = 0
        skipping  = False # skip first time step of every file except the very first file (i.e. ignore repetitive data)
        with open(file_name+"_"+str(i)+".dat", "r") as f:
            for line in f:
                if len(line.strip()) != 0 : # make sure line isn't empty
                    line_split = line.split()
                    if line[0] in ['#']: # title, labels, and comments in file
                        if line_split[1] == 'Time':
                            if n_skipped == 0 and i != 1:
                                skipping = True
                                n_skipped += 1
                            else:
                                time_step       = int(float(line_split[3]))
                                data[time_step] = copy.deepcopy(dict_temp)
                                skipping        = False
                        continue
                    elif not skipping:
                        for j in range(len(parameters)):
                            value = float(line_split[j])
                            if parameters[j] == "twist":
                                value = normalize_range(value, 0, 360)
                            data[time_step][parameters[j]].append(value)
    return data

################################################################################################
#
# MAIN PROGRAM
#
################################################################################################

def main():
    data_per_dt = 20
    data = get_data(paths[0]+"L-BPS")
    print(data[0]["twist"])
    print(data[30000]["twist"])
    print(data[30050]["twist"])
    #for x in data:
    #    print(x)
    #    for y in data[x]:
    #        print(y)
            #print(y,':',data[x][y])

if __name__ == "__main__": 
    main()
            