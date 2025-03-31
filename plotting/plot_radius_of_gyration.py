"""
    Author: Rachel Bricker
    Year:   2025
"""

"""
    Plots the radius of gyration as a function of time.
"""

"""
   usage: python3 plot_radius_of_gyration.py
      i. [path to .xvg file outputted by GROMACS utility `gyrate`]
   
   example: python3 plot_radius_of_gyration.py \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/ssDNA1/gyrate.xvg \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/ssDNA2/gyrate.xvg \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/ssDNA3/gyrate.xvg \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/ssDNA4/gyrate.xvg \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/ssDNA5/gyrate.xvg
"""

import sys
import numpy as np 
from functions_for_plots import *

# command line input
paths = list(sys.argv[1:])

################################################################################################
#
# FUNCTIONS
#
################################################################################################

def get_time_and_gyrate(paths):
    gyrate = [[] for i in range(len(paths))]
    time   = []

    # iterate over .xvg files passed through command line
    for i in range(len(paths)):
        file = paths[i]
        with open(file, "r") as f:
            data = read_xvg_file(f) # read data file

            # iterate over configurations
            for j in range(len(data)):
                if i == 0:
                    time.append(data[j][0]) # get time (measured in ps)

                gyrate[i].append(data[j][1]) # get radius of gyration of DNA

    time = [t/1000 for t in time] # convert ps -> ns

    return time, gyrate

################################################################################################
#
# MAIN PROGRAM
#
################################################################################################

def main():
    time, gyrate = get_time_and_gyrate(paths)
    font_leg     = set_rcParameters()

    # set figure dimensions
    fig_height = 1.55   
    fig_width  = 2.825
    fig, ax    = plt.subplots(1, figsize=(fig_width, fig_height))

    # plot data
    for i in range(len(gyrate)):
        plt.plot(time, moving_average(moving_average(gyrate[i], 500), 100))

    # set x-axis label
    plt.xlabel("Time, ns")
        
    # set y-axis label
    plt.ylabel("Radius of Gyration, nm")

    # set x-axis limits
    plt.xlim(-8, 308)

    # show grid
    plt.grid()

    plt.tight_layout()

    # save figure
    plt.savefig("fig4_gyrate.svg", bbox_inches="tight", dpi=600)

if __name__ == "__main__": 
    main()
    