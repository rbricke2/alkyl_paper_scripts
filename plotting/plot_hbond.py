"""
    Author: Rachel Bricker
    Year:   2025
"""

"""
    Plots hydrogen bond analysis.
"""

"""
   usage: python3 plot_radius_of_gyration.py
      1. [list of scenario names for legend (must be parallel w.r.t the arguments i, ..., i+n)]
      2. [name of .xvg file outputted by GROMACS utility `distance`]
      3. [name of .xvg file outputted by GROMACS utility `angle`]
      i. [path to directories that contain arguments (2.) and (3.) that you want to plot]
   
   example: python3 plot_hbond.py \
            "unmodified,4 decyls\n(together),4 decyls\n(spaced out),10 decyls,10 ethyls" \
            hbond.xvg \
            hbond_angle.xvg \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/AMBER/dsDNA1/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/AMBER/dsDNA2/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/AMBER/dsDNA3/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/AMBER/dsDNA4/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/AMBER/dsDNA5/
"""

import sys
import numpy as np 
from functions_for_plots import *
import matplotlib as mpl

# command line input
input_list = sys.argv[1].replace("\\n", "\n").replace("\\(", "(").replace("\\)", ")")
legend     = input_list.split(',')
dist_xvg   = str(sys.argv[2])
ang_xvg    = str(sys.argv[3])
paths      = list(sys.argv[4:])

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

def plot_color_map(time, hbond_bool_matrix, font_leg):
    font_size   = font_leg.get_size()
    font_family = font_leg.get_family()[0]

    # figure dimensions  
    fig_width  = 8.7
    fig_height = (fig_width/7.5)

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

        # set tick frequency on x- and y-axis
        axes[scenario].set_yticks(np.arange(ybottom, ytop, 100))
        axes[scenario].set_xticks(np.arange(xleft, xright, 4))

        if scenario == len(hbond_bool_matrix)//2:
            axes[scenario].set_xlabel("Base pair")

    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)

    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)

    # set y-axis label
    plt.ylabel("Time (ns)")

    # color bar
    cbar = fig.colorbar(pixel_plot, ax=axes.ravel().tolist(), ticks=[0,1], pad=0.02, aspect=10)
    cbar.ax.set_yticklabels(['Hydrogen\nbonding', 'Broken\nhydrogen\nbonding'])

    # change font size for color bar
    cbar.ax.tick_params(labelsize=font_size)

    # change font for color bar
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_family(font_family)

    #plt.tight_layout()

    # save figure
    plt.savefig("colormap_hbond_binary.svg", bbox_inches="tight", dpi=600)

def plot_data(x, y, x_label, y_label, title, legend, file_name, fig_width, fig_height):
    # set rcParams
    font_leg = set_rcParameters()

    # set figure dimensions
    fig, ax     = plt.subplots(1, figsize=(fig_width, fig_height))

    # plot data
    for i in range(len(y)):
        plt.plot(x, moving_average(moving_average(y[i], 500), 100), label=legend[i])

    # position legend to the left
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop=font_leg) 

    # set x-axis label
    plt.xlabel(x_label)
        
    # set y-axis label
    plt.ylabel(y_label)

    # set title
    if title != None:
        plt.title(title)

    # show grid
    plt.grid()

    plt.tight_layout()

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
    plot_color_map(time, hbond_bool_matrix, font_leg)
    
    # plot other data as function of time
    x_label     = "Time, ns"
    fig_width   = 5
    golden_mean = (np.sqrt(5)-1.0)/2.0     # aesthetic ratio
    fig_height  = fig_width*golden_mean    # height in inches
    plot_data(time,
              n_broken_hbond,
              x_label,
              "Number of broken\nhydrogen bonds",
              "Entire Duplex (Parmbsc1)",
              legend,
              "melted_hbond_vs_time.svg",
              fig_width,
              fig_height)
    plot_data(time,
              get_avg_dist_per_conf(distances, stop_residue_id),
              x_label,
              r"$\langle r \rangle$, nm",
              "First Six Base Pairs (Parmbsc1)",
              legend,
              "distance_vs_time_terminal.svg",
              fig_width,
              fig_height)
    plot_data(time,
              get_n_broken_hbond(hbond_bool_matrix, stop_residue_id),
              x_label,
              "Number of broken\nhydrogen bonds",
              "First Six Base Pairs (Parmbsc1)",
              legend,
              "melted_hbond_vs_time_terminal.svg",
              fig_width,
              fig_height)
    
    # print statistics
    for i in range(len(n_broken_hbond)):
        print("Average number of melted hydrogen bonds for file " + str(i+1) + " (excluding first 200 ns): " + str(round(statistics.mean(n_broken_hbond[i][4000:]),3)) + " +/- " + str(round(statistics.stdev(n_broken_hbond[i][4000:]),3)))

if __name__ == "__main__": 
    main()
    