"""
    Author: Rachel Bricker
    Year:   2025
"""

"""
    Plots hydrogen bond analysis.
"""

"""
   usage: python3 plot_hbond.py
      1. list of model names for legend (must be parallel w.r.t the arguments i to i+n)
      2. name of .xvg file outputted by GROMACS utility `distance`
      3. name of .xvg file outputted by GROMACS utility `angle`
      4. enter 0 if constant temperature, 1 if simulated annealing
      i. path to directories that contain arguments (2.) and (3.) that you want to plot
   
   examples: python3 plot_hbond.py \
             "(a),(b),(d)" \
             hbond.xvg \
             hbond_angle.xvg \
             0 \
             /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/AMBER/dsDNA1/ \
             /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/AMBER/dsDNA3/ \
             /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/AMBER/dsDNA4/
             
             python3 plot_hbond.py \
             "(a),(b),(d)" \
             hbond.xvg \
             hbond_angle.xvg \
             1 \
             /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/annealing_AMBER/dsDNA1/ \
             /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/annealing_AMBER/dsDNA3/ \
             /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/annealing_AMBER/dsDNA4/
"""

import sys
import numpy as np 
from functions_for_plots import *
import matplotlib as mpl
from matplotlib.ticker import FixedLocator

# command line input
input_list = sys.argv[1].replace("\\n", "\n").replace("\\t", "\t").replace("\\(", "(").replace("\\)", ")")
legend     = input_list.split(',')
dist_xvg   = str(sys.argv[2])
ang_xvg    = str(sys.argv[3])
annealing  = bool(int(sys.argv[4]))
paths      = list(sys.argv[5:])

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

def get_dist(paths, dist_xvg):
    time      = []
    distances = [ [] for scenario in range(len(paths)) ]
    
    # loop over scenarios
    for scenario in range(len(paths)):
        file = paths[scenario] + dist_xvg
        with open(file, "r") as f:
            data = read_xvg_file(f) # read data file

            # loop over configurations
            for conf in range(len(data)):
                time_step = data[conf].pop(0)
                
                if scenario == 0:
                    time.append(time_step) # get time (measured in ps)

                # record distances
                distances[scenario].append(data[conf])
    
    time = [t/1000 for t in time] # convert ps -> ns
    
    return time, distances

def get_angle(paths, ang_xvg):
    angles = [ [] for scenario in range(len(paths)) ]
    
    # loop over scenarios
    for scenario in range(len(paths)):
        file = paths[scenario] + ang_xvg
        with open(file, "r") as f:
            data = read_xvg_file(f) # read data file

            # loop over configurations
            for conf in range(len(data)):
                data[conf].pop(0)

                # record distances
                angles[scenario].append(data[conf])
    
    return angles

def get_hbond_existence(distances, angles):
    # hbond exists if:
    #     distance <= 0.35 nm
    #     angle    <= 30 degrees

    # list of matrices
    hbond_bool_matrix = []

    for scenario in range(len(distances)):
        # initialize data matrix
        Z = []

        # number of configurations
        n_configurations = len(distances[scenario])

        # get number of base pairs
        base_pairs   = list(range(1, len(distances[scenario][0])+1))
        n_base_pairs = len(base_pairs)
        del base_pairs

        # loop over configurations
        for conf in range(n_configurations):
            if not len(Z): # if data matrix is empty...
                Z = np.zeros((n_configurations, n_base_pairs)) # ... then create matrix of zeros

            # loop over base pairs
            for base_pair in range(n_base_pairs):
                dist  = distances[scenario][conf][base_pair]
                angle = angles[scenario][conf][base_pair]

                if not (dist <= 0.35 and angle <= 30):
                    Z[conf][base_pair]  = 1
        
        # append matrix
        hbond_bool_matrix.append(Z)

    return hbond_bool_matrix

def get_n_broken_hbond(hbond_bool_matrix, stop_residue_id=None):
    # hbond exists if:
    #     distance <= 0.35 nm
    #     angle    <= 30 degrees
    
    n_broken_hbond = [ [] for scenario in range(len(hbond_bool_matrix)) ]
    
    # loop over scenarios
    for scenario in range(len(hbond_bool_matrix)):

        # loop over configurations
        for conf in range(len(hbond_bool_matrix[scenario])):
            broken_hbond_count = 0
            
            if stop_residue_id == None:
                stop_residue_id = len(hbond_bool_matrix[scenario][conf])

            # loop over base pairs
            for base_pair in range(stop_residue_id):
                if hbond_bool_matrix[scenario][conf][base_pair] == 1:
                    broken_hbond_count += 1

            n_broken_hbond[scenario].append(broken_hbond_count)

    return n_broken_hbond
    
def get_avg_dist_per_conf(distances, stop_residue_id=None):  
    dist_avg = [ [] for scenario in range(len(distances)) ]
    
    # loop over scenarios
    for scenario in range(len(distances)):

        # loop over configurations
        for conf in range(len(distances[scenario])):
            
            if stop_residue_id == None:
                stop_residue_id = len(distances[scenario][conf])

            dist_avg[scenario].append(statistics.mean(distances[scenario][conf][:stop_residue_id]))
    
    return dist_avg

def plot_color_map(time, hbond_bool_matrix, annealing, font_leg):
    font_size   = font_leg.get_size()
    font_family = font_leg.get_family()[0]
    
    # initialize file name
    file_name = "colorplot_hbond_binary.svg"

    # figure dimensions  
    fig_width  = 8.7
    fig_height = (fig_width/7.5)
    fig_height = 2.2
    fig_width  = (len(hbond_bool_matrix)/5)*(fig_width)
    fig_width  = 4.1
    
    if annealing:
        fig_width = 4.1
        time      = [t+600 for t in time]
    
    # initialize color bar padding
    padding = 0.03

    # initialize figure
    fig, axes  = plt.subplots(nrows=1, ncols=len(hbond_bool_matrix), sharey=True, figsize=(fig_width, fig_height))

    # helpful for 2d color plots: https://www.geeksforgeeks.org/create-2d-pixel-plot-in-python/
    #cmap = mpl.colors.ListedColormap([(0.827, 0.933, 1), (0.827, 0.161, 0.102)])
    cmap = mpl.colors.ListedColormap([(0.922, 0.922, 0.922), (0.62, 0.192, 0.961)])

    for scenario in range(len(hbond_bool_matrix)):
        Z = hbond_bool_matrix[scenario]

        time_step = time[1]-time[0]
        xleft     = 1
        xright    = np.shape(Z)[1]+1
        ybottom   = time[0]
        ytop      = time[-1]+time_step

        #ratio         = 1.0
        #aspect_square = abs((xright-xleft)/(ybottom-ytop))*ratio

        # plot color plot
        pixel_plot = axes[scenario].imshow(Z, cmap=cmap, interpolation='nearest', origin='lower', vmin=0, vmax=1,
                                           extent=[xleft, xright, ybottom, ytop], aspect='auto')

        axes[scenario].set_title(legend[scenario])

        # set tick frequency on x- and y-axis and center ticks on row/column
        half_time_step = time_step/2
        x_axis_freq = 5
        y_axis_freq = 100
        if annealing:
            y_axis_freq = 200
        axes[scenario].set_yticks(np.arange(ybottom+half_time_step, ytop+half_time_step, y_axis_freq), np.arange(ybottom, ytop, y_axis_freq, dtype=int))
        axes[scenario].set_xticks(np.arange(xleft+0.5, xright+0.5, x_axis_freq), np.arange(xleft, xright, x_axis_freq))
        
        # set minor tick locations on the x-axis
        axes[scenario].xaxis.set_minor_locator(FixedLocator(np.arange(xleft+0.5, xright+0.5, 1)))

        # put x-axis label on centermost plot
        if scenario == len(hbond_bool_matrix)//2:
            axes[scenario].set_xlabel("Base pair")
        # put y-axis label on rightmost plot
        if scenario == 0:            
            axes[scenario].set_ylabel("Simulation time (ns)")

        if annealing: 
            # add secondary y-axis
            ax2 = axes[scenario].secondary_yaxis('right', functions=(lambda x: 300+((x-600)*((440-300)/1200)), lambda x: 300+((x-600)*((440-300)/1200))))
            
            # set tick frequency on secondary y-axis and center ticks on row/column
            freq = 20
            ax2.set_yticks(np.arange(300, 440+freq, freq, dtype=int))
            
            # add extra padding to color bar to create room for secondary y-axis
            padding += 0.033
            
            if scenario == (len(hbond_bool_matrix)-1):            
                # label secondary y-axis on leftmost plot
                ax2.set_ylabel('Temperature (K)', rotation=270, va='bottom')
            else:
                # remove secondary y-axis tick labels on all plots except for leftmost
                ax2.set_yticklabels([])
            
            # change file name
            file_name = "colorplot_hbond_binary_annealing.svg"
    

    # color bar
    cbar = fig.colorbar(pixel_plot, ax=axes, ticks=[0,1], pad=0.3, fraction=0.05, aspect=30, orientation='horizontal')
    cbar.ax.set_xticklabels(['Hydrogen bonding', 'Broken hydrogen bonding'])

    # change font size for color bar
    cbar.ax.tick_params(labelsize=font_size)

    # change font for color bar
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_family(font_family)

    #plt.tight_layout()

    # save figure
    plt.savefig(file_name, bbox_inches="tight", dpi=600)

################################################################################################
#
# MAIN PROGRAM
#
################################################################################################

def main():
    stop_residue_id = 6   # first six residues
    
    # get data from .xvg files
    time, distances   = get_dist(paths, dist_xvg)
    angles            = get_angle(paths, ang_xvg)
    hbond_bool_matrix = get_hbond_existence(distances, angles)
    n_broken_hbond    = get_n_broken_hbond(hbond_bool_matrix)

    # set rcParams
    font_leg = set_rcParameters()

    # plot boolean color map
    plot_color_map(time, hbond_bool_matrix, annealing, font_leg)

    """
    # plot other data as function of time
    x_label     = "Simulation time (ns)"
    fig_width   = 5
    golden_mean = (np.sqrt(5)-1.0)/2.0     # aesthetic ratio
    fig_height  = fig_width*golden_mean    # height in inches
    plot_data(time,
              n_broken_hbond,
              x_label,
              "Number of broken base pairs (bp)",
              "Entire Duplex (Parmbsc1)",
              legend,
              "melted_hbond_vs_time.svg",
              fig_width,
              fig_height)
    plot_data(time,
              get_avg_dist_per_conf(distances, stop_residue_id),
              x_label,
              r"$\langle r \rangle$ (nm)",
              "First Six Base Pairs (Parmbsc1)",
              legend,
              "distance_vs_time_terminal.svg",
              fig_width,
              fig_height)
    plot_data(time,
              get_n_broken_hbond(hbond_bool_matrix, stop_residue_id),
              x_label,
              "Number of broken base pairs (bp)",
              "First Six Base Pairs (Parmbsc1)",
              legend,
              "melted_hbond_vs_time_terminal.svg",
              fig_width,
              fig_height)
    """

    # print statistics
    frame_200ns = 4000
    frame_600ns = 12000
    for i in range(len(n_broken_hbond)):
        print("Average number of melted base pairs for file " + str(i+1) + ": " + str(round(statistics.mean(n_broken_hbond[i]),1)) + " +/- " + str(round(statistics.stdev(n_broken_hbond[i]),1)))
    for i in range(len(n_broken_hbond)):
        if len(n_broken_hbond[i]) > frame_200ns+1:
            print("Average number of melted base pairs for file " + str(i+1) + " (excluding first 200 ns): " + str(round(statistics.mean(n_broken_hbond[i][frame_200ns:]),1)) + " +/- " + str(round(statistics.stdev(n_broken_hbond[i][frame_200ns:]),1)))
    for i in range(len(n_broken_hbond)):
        if len(n_broken_hbond[i]) > frame_600ns+1:
            print("Average number of melted base pairs for file " + str(i+1) + " (excluding first 600 ns): " + str(round(statistics.mean(n_broken_hbond[i][frame_600ns:]),1)) + " +/- " + str(round(statistics.stdev(n_broken_hbond[i][frame_600ns:]),1)))

if __name__ == "__main__": 
    main()
    