"""
    Author: Rachel Bricker
    Year:   2025
"""

"""
    Plots the radius of gyration as a function of time.
"""

"""
   usage: python3 plot_radius_of_gyration.py
      1. figure width
      2. figure height
      i. path to .xvg file outputted by GROMACS utility `gyrate` that you want to plot
   
   examples: python3 plot_radius_of_gyration.py \
             3.7 \
             1.57 \
             /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/ssDNA1/gyrate.xvg \
             /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/ssDNA3/gyrate.xvg \
             /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/ssDNA2/gyrate.xvg \
             /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/ssDNA4/gyrate.xvg \
             /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/CHARMM36/ssDNA5/gyrate.xvg
             
             python3 plot_radius_of_gyration.py \
             7.1 \
             2.5 \
             /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/AMBER/ssDNA1/gyrate.xvg \
             /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/AMBER/ssDNA3/gyrate.xvg \
             /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/AMBER/ssDNA2/gyrate.xvg \
             /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/AMBER/ssDNA4/gyrate.xvg \
             /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/AMBER/ssDNA5/gyrate.xvg
"""

import sys
import numpy as np 
from functions_for_plots import *

# command line input
fig_width  = float(sys.argv[1])
fig_height = float(sys.argv[2])
paths      = list(sys.argv[3:])

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
    fig, ax = plt.subplots(1, figsize=(fig_width, fig_height))

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
    plt.savefig("gyrate_plot.svg", bbox_inches="tight", dpi=600)
    
    # print statistics
    frame_200ns  = 4000
    rounding     = 2
    for i in range(len(gyrate)):
        print("Average radius of gyration value for file " + str(i+1) + ": " + str(round(statistics.mean(gyrate[i]),rounding)) + " +/- " + str(round(statistics.stdev(gyrate[i]),rounding)))
        
    for i in range(len(gyrate)):
        print("Average radius of gyration value for file " + str(i+1) + " (excluding first 200 ns): " + str(round(statistics.mean(gyrate[i][frame_200ns:]),rounding)) + " +/- " + str(round(statistics.stdev(gyrate[i][frame_200ns:]),rounding)))

if __name__ == "__main__": 
    main()
    