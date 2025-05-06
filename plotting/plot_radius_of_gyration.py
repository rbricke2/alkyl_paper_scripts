"""
    Author: Rachel Bricker
    Year:   2025
"""

"""
    Plots the radius of gyration as a function of time.
"""

"""
   usage: python3 plot_radius_of_gyration.py
      i. path to .xvg file outputted by GROMACS utility `gyrate` that you want to plot
   
   example: python3 plot_radius_of_gyration.py \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/ssDNA1/gyrate.xvg \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/ssDNA3/gyrate.xvg \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/ssDNA2/gyrate.xvg \
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

    # iterate over each .xvg file passed via command line
    for i in range(len(paths)):
        file = paths[i]
        with open(file, "r") as f:
            data = read_xvg_file(f) # read data file

            # iterate over configurations
            for j in range(len(data)):
                if i == 0:
                    time.append(data[j][0]) # get time (measured in ps)

                gyrate[i].append(data[j][1]) # get radius of gyration of molecule

    time = [t/1000 for t in time] # convert ps -> ns

    return time, gyrate

################################################################################################
#
# MAIN PROGRAM
#
################################################################################################

def main():
    # get data from .xvg files
    time, gyrate = get_time_and_gyrate(paths)
    
    # set rcParams
    font_leg = set_rcParameters()

    # set figure dimensions
    fig_height = 1.57   
    fig_width  = 3
    fig, ax    = plt.subplots(1, figsize=(fig_width, fig_height))

    # plot data
    for i in range(len(gyrate)):
        plt.plot(time, moving_average(moving_average(gyrate[i], 500), 100))

    # set x-axis label
    plt.xlabel("Simulation time (ns)")
        
    # set y-axis label
    plt.ylabel("Gyration radius (nm)")

    # set x-axis limits
    plt.xlim(time[0]-8, time[-1]+8)

    # show grid
    plt.grid()

    plt.tight_layout()

    # save figure
    plt.savefig("fig4_gyrate.svg", bbox_inches="tight", dpi=600)
    
    # print statistics
    for i in range(len(gyrate)):
        print("Average radius of gyration value for file " + str(i+1) + ": " + str(round(statistics.mean(gyrate[i]),3)) + " +/- " + str(round(statistics.stdev(gyrate[i]),3)))
        
    for i in range(len(gyrate)):
        print("Average radius of gyration value for file " + str(i+1) + " (excluding first 200 ns): " + str(round(statistics.mean(gyrate[i][4000:]),3)) + " +/- " + str(round(statistics.stdev(gyrate[i][4000:]),3)))

if __name__ == "__main__": 
    main()
    