"""
    Author:          Rachel Bricker
    Year:            2025
"""

"""
   usage: python3 plot_x3DNA.py
      1. list of model names for legend (must be parallel w.r.t the arguments i to i+n)
      2. file name
      3. parameter
      i. path to directories that contain arguments (2.) and (3.) that you want to plot

   examples: python3 plot_x3DNA.py \
            "(a),(b),(c),(d),(e)" \
            L-BPS \
            twist \
            "Average twist, degrees" \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/movies_gifs_pictures/dsDNA1_CHARMM/x3DNA/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/movies_gifs_pictures/dsDNA3_CHARMM/x3DNA/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/movies_gifs_pictures/dsDNA2_CHARMM/x3DNA/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/movies_gifs_pictures/dsDNA4_CHARMM/x3DNA/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/movies_gifs_pictures/dsDNA5_CHARMM/x3DNA/
            
            python3 plot_x3DNA.py \
            "(a),(b),(c),(d),(e)" \
            L-BPS \
            twist \
            "Average twist, degrees" \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/movies_gifs_pictures/dsDNA1/x3DNA/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/movies_gifs_pictures/dsDNA3/x3DNA/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/movies_gifs_pictures/dsDNA2/x3DNA/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/movies_gifs_pictures/dsDNA4/x3DNA/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/movies_gifs_pictures/dsDNA5/x3DNA/
"""

import sys
import numpy as np 
import math
import copy
from functions_for_plots import *

# command line input
input_list  = sys.argv[1].replace("\\n", "\n").replace("\\t", "\t").replace("\\(", "(").replace("\\)", ")")
legend      = input_list.split(',')
file_prefix = str(sys.argv[2])
parameter   = str(sys.argv[3])
y_label     = str(sys.argv[4])
paths       = list(sys.argv[5:])

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
    dict_temp = {}

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
                        else:
                            dict_temp = {str(line_split[0][1:]).lower(): [],
                                         str(line_split[1]).lower(): [],
                                         str(line_split[2]).lower(): [],
                                         str(line_split[3]).lower(): [],
                                         str(line_split[4]).lower(): [],
                                         str(line_split[5]).lower(): []}
                            # list of DNA parameters
                            parameters = list(dict_temp.keys())
                        continue
                    elif not skipping:
                        for j in range(len(parameters)):
                            param = parameters[j]
                            value = float(line_split[j])
                            #if param == "twist":
                                #value = normalize_range(value, 0, 360) # convert range of angle from [-180, 180] -> [0, 360]
                            data[time_step][param].append(value)
    return data

# this function was grabbed from:
# https://matplotlib.org/stable/gallery/subplots_axes_and_figures/secondary_axis.html
def forward_inverse(x):
    x = np.array(x, float)
    near_zero = np.isclose(x, 0)
    x[near_zero] = np.inf
    x[~near_zero] = 360 / x[~near_zero]
    return x

################################################################################################
#
# MAIN PROGRAM
#
################################################################################################

def main():
    # set rcParams
    font_leg = set_rcParameters()

    # set figure dimensions
    fig, ax1 = plt.subplots(1, figsize=(7, 1.8))

    for i in range(len(paths)):
        data = get_data(paths[i]+file_prefix)

        # get the average twist per configuration
        avg_param = []
        for dt in data:
            if parameter == "twist":
                # exclude terminal base pairs at each end of the DNA duplex for twist
                n_bp_excluded = 3
                avg_param.append(statistics.mean(data[dt][parameter][n_bp_excluded:-n_bp_excluded]))
            else:
                avg_param.append(statistics.mean(data[dt][parameter]))

        # plot data
        ax1.plot(list(data.keys()), moving_average(moving_average(avg_param, 500), 100), label=legend[i])

    # position legend to the left
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop=font_leg)
    plt.legend(loc='best')
    
    # set x-axis label
    plt.xlabel("Simulation time, ns")

    # label right y-axis for the middle row
    ax1.set_ylabel(y_label)
    
    if parameter == "twist":
    # secondary axis
        ax2 = ax1.secondary_yaxis('right', functions=(forward_inverse, forward_inverse))

        # label left y-axis for the middle row
        ax2.set_ylabel('Helical repeat, bp/turn', rotation=270, va='bottom')

    # show grid
    plt.grid()

    plt.tight_layout()

    # save figure
    plt.savefig(parameter+"_plot.svg", bbox_inches="tight", dpi=600)

if __name__ == "__main__": 
    main()
            