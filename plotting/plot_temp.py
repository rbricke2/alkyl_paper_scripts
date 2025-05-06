"""
    Author: Rachel Bricker
    Year:   2025
"""

"""
    Plots temperature.
"""

"""
   usage: python3 plot_temp.py
      1. list of scenario names for legend (must be parallel w.r.t the arguments i to i+n)
      i. path to .xvg file outputted by GROMACS utility `energy` that you want to plot
   
   example: python3 plot_temp.py \
            "(a),(b),(d)" \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/annealing_AMBER/dsDNA1/temperature_1200ns.xvg \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/annealing_AMBER/dsDNA3/temperature_1200ns.xvg \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/annealing_AMBER/dsDNA4/temperature_1200ns.xvg 
"""

import sys
import numpy as np 
from functions_for_plots import *
import matplotlib as mpl
from matplotlib.ticker import FixedLocator

# command line input
input_list = sys.argv[1].replace("\\n", "\n").replace("\\(", "(").replace("\\)", ")")
legend     = input_list.split(',')
paths      = list(sys.argv[2:])

################################################################################################
#
# FUNCTIONS
#
################################################################################################

def get_time_and_temp(paths):
    temp = [[] for i in range(len(paths))]
    time = []

    # iterate over each .xvg file passed via command line
    for i in range(len(paths)):
        file = paths[i]
        with open(file, "r") as f:
            data = read_xvg_file(f) # read data file

            # iterate over configurations
            for j in range(len(data)):
                if i == 0:
                    time.append(data[j][0]) # get time (measured in ps)

                temp[i].append(data[j][1]) # get temperature of system (measured in K)

    time = [t/1000 for t in time] # convert ps -> ns

    return time, temp

################################################################################################
#
# MAIN PROGRAM
#
################################################################################################

def main():
    # get data from .xvg files
    time, temp = get_time_and_temp(paths)
    
    # set rcParams
    font_leg = set_rcParameters()
    
    # set figure dimensions
    fig_width   = 5
    golden_mean = (np.sqrt(5)-1.0)/2.0     # aesthetic ratio
    fig_height  = fig_width*golden_mean    # height in inches

    # plot the data
    plot_data(time,
              temp,
              "Simulation time (ns)",
              "Temperature (K)",
              "Simulated Annealing",
              legend,
              "temp_vs_time.svg",
              fig_width,
              fig_height)

if __name__ == "__main__": 
    main()
    