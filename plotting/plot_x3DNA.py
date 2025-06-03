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
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/movies_gifs_pictures/dsDNA1_CHARMM/x3DNA/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/movies_gifs_pictures/dsDNA3_CHARMM/x3DNA/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/movies_gifs_pictures/dsDNA2_CHARMM/x3DNA/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/movies_gifs_pictures/dsDNA4_CHARMM/x3DNA/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/movies_gifs_pictures/dsDNA5_CHARMM/x3DNA/
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
    # template dictionary
    dict_temp = {"shift": [], "slide": [], "rise": [], "tilt": [], "roll": [], "twist": []}

    # list of DNA parameters
    parameters = list(dict_temp.keys())

    # stores the data
    data = {}

    # records the time step
    time_step = None

    # iterating over 20 files because I analyzed the 600 ns trajectory in 30 ns intervals
    for i in range(1,21):
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
                                time_step       = float(line_split[3])/1000 # convert ps -> ns
                                data[time_step] = copy.deepcopy(dict_temp)
                                skipping        = False
                        continue
                    elif not skipping:
                        for j in range(len(parameters)):
                            param = parameters[j]
                            value = float(line_split[j])
                            if param == "twist":
                                value = normalize_range(value, 0, 360) # convert range of angle from [-180, 180] -> [0, 360]
                            data[time_step][param].append(value)
    return data

################################################################################################
#
# MAIN PROGRAM
#
################################################################################################

def main():
    # set rcParams
    font_leg = set_rcParameters()

    # set figure dimensions
    fig, ax = plt.subplots(1, figsize=(8.2, 5))

    for i in range(len(paths)):
        data = get_data(paths[i]+"L-BPS")

        # get the average twist per configuration
        avg_twist = []
        for dt in data:
            # exclude first two terminal base pairs
            n_exclude = 2
            avg_twist.append(statistics.mean(data[dt]["twist"][n_exclude:-n_exclude]))

        # plot data
        plt.plot(list(data.keys()), moving_average(moving_average(avg_twist, 500), 100), label=legend[i])

    # position legend to the left
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop=font_leg) 

    # set x-axis label
    plt.xlabel("Simulation time, ns")

    # set y-axis label
    plt.ylabel("Average twist, degrees")

    # show grid
    plt.grid()

    plt.tight_layout()

    # save figure
    plt.savefig("twist_plot.png", bbox_inches="tight", dpi=600)

if __name__ == "__main__": 
    main()
            