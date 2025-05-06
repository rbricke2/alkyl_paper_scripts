"""
    Author: Rachel Bricker
    Year:   2025
"""

"""
    Plots stacking analysis.
"""

"""
   usage: python3 plot_stacking.py
      1. 
   
   example: python3 plot_stacking.py \
            "(a),(b),(d)" \
            com_files \
            vec_files \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/annealing_AMBER/dsDNA1/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/annealing_AMBER/dsDNA3/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/annealing_AMBER/dsDNA4/ 
"""

import sys
import numpy as np 
from functions_for_plots import *

# command line input
input_list = sys.argv[1].replace("\\n", "\n").replace("\\(", "(").replace("\\)", ")")
legend     = input_list.split(',')
com_dir    = str(sys.argv[2])
vec_dir    = str(sys.argv[3])
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


def S(alpha):
    return ( np.exp( -(alpha**4) ) +
             np.exp( -((alpha - np.pi)**4) ) +
             (0.1*np.exp( -((alpha - (0.5*np.pi))**4) )) )

# stacking coordinate
def xi(COM_dist, alpha):
    return ( COM_dist / ( S(alpha) ) )

def get_data(paths, com_dir, vec_dir):
    stacking_coords = [ [] for scenario in range(len(paths)) ]
    for scenario in range(len(paths)):
        # check if the DNA is double-stranded or single-stranded
        ds = False
        if "ds" in paths[scenario].split("/")[-2]:
            ds         = True
            n_residues = 42
        else:
            n_residues = 21

        time       = []
        COM_coords = [ [] for resi in range(n_residues) ]

        # iterate over residues
        for resi in range(n_residues):
            file = paths[scenario] + com_dir + "/nucleobase_COM_coord_" + str(resi+1) + ".xvg"
            with open(file, "r") as f:
                data = read_xvg_file(f) # read data file

                # iterate over time
                for t in range(len(data)):
                    time_step = data[t].pop(0)

                    # only record time once
                    if resi == 0:
                        time.append(time_step) # get time (measured in ps)

                    COM_coords[resi].append(data[t])

        time = [t/1000 for t in time] # convert ps -> ns

        atoms_xyz = []
        norms     = [ [] for resi in range(n_residues) ] # base plane normal vectors

        # iterate over residues
        for resi in range(n_residues):
            file = paths[scenario] + vec_dir + "/nucleobase_vec_coord_" + str(resi+1) + ".xvg"
            with open(file, "r") as f:
                data = read_xvg_file(f) # read data file

                # iterate over time
                for t in range(len(data)):
                    data[t].pop(0) # time was already recorded, so discard it

                    # record the xyz coordinates of atoms
                    # data is written like: a1x, a1y, a1z, a2x, a2y, a2z, ...
                    atoms_xyz = [data[t][i:i+3] for i in range(0, len(data[t]), 3)]

                    COM_coord = COM_coords[resi][t]
                    vec_a     = np.array([atoms_xyz[0][0]-COM_coord[0],
                                          atoms_xyz[0][1]-COM_coord[1],
                                          atoms_xyz[0][2]-COM_coord[2]])
                    vec_b     = np.array([atoms_xyz[1][0]-COM_coord[0],
                                          atoms_xyz[1][1]-COM_coord[1],
                                          atoms_xyz[1][2]-COM_coord[2]])
                                          
                    # record normal vector of base plane
                    norms[resi].append(np.cross(vec_a, vec_b))

        # measure the angle alpha, in radians, between base planes
        # and the distance between mass centers of consecutive bases

        # iterate over time
        for t in range(len(time)):
            time_step = []
            # iterate over residues
            for resi in range(n_residues):
                if ( (ds and resi == (n_residues//2)) or (resi+1 == n_residues) ):
                    # we are at a terminal, so there is no following base
                    continue

                # distance between mass centers
                first_COM_pos  = COM_coords[resi][t]      # grab the center of mass position of nucleobase i
                second_COM_pos = COM_coords[resi+1][t]    # grab the center of mass position of nucleobase i+1

                # get the differences between the corresponding x, y, and z values for the distance formula
                delta_x = second_COM_pos[0] - first_COM_pos[0]    # change in x-coordinate
                delta_y = second_COM_pos[1] - first_COM_pos[1]    # change in y-coordinate
                delta_z = second_COM_pos[2] - first_COM_pos[2]    # change in z-coordinate

                # apply distance formula
                COM_dist = np.sqrt(delta_x**2 + delta_y**2 + delta_z**2)

                # angle between base planes
                u         = norms[resi][t]      # normal of base plane of residue i
                v         = norms[resi+1][t]    # normal of base plane of residue i+1
                cos_theta = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
                alpha     = np.arccos(np.clip(cos_theta, -1.0, 1.0))    # angle between consecutive base planes (radians)

                time_step.append(xi(COM_dist, alpha))
            
            stacking_coords[scenario].append(time_step)    # the stacking coordinate, xi, is measured in nm

    return time, stacking_coords

def analyze_data(stacking_coords, max_strand_len=21):
    # records the number of consecutively stacked bases for each scenario for each configuration
    consecutive_stacked = [ [ [] for t in range(len(stacking_coords[scenario])) ] for scenario in range(len(stacking_coords)) ]
    
    # records the number of broken stacked bases for each scenario for each configuration
    n_broken_stacking = [ [] for scenario in range(len(stacking_coords)) ]
    
    # will allow us to determine if we are at the last recorded stacking coordinate of a strand
    n_stacking_coord_per_strand = max_strand_len-1
    
    # transient point, measured in nm, that separates the pools of the stacked and unstacked configurations
    transient_pt = 0.6

    # loop over scenarios
    for scenario in range(len(stacking_coords)):
        # loop over time
        for t in range(len(stacking_coords[scenario])):

            # iterate through the stacking coordinate values
            consecutive_count = 1   # tracks the number of consecutive bases stacked
            broken_count      = 0   # tracks the number of broken stacking
            prev              = -1  # records previous data point

            # loop over recorded stacking coordinates at time 't'
            for coord_itr in range(len(stacking_coords[scenario][t])):
                coord_value = stacking_coords[scenario][t][coord_itr]
                if coord_value > transient_pt:
                    if prev <= transient_pt and prev != -1:
                        consecutive_stacked[scenario][t].append(consecutive_count)
                    consecutive_count = 1
                    broken_count     += 1
                elif (coord_itr == n_stacking_coord_per_strand-1) or (coord_itr == (n_stacking_coord_per_strand*2)-1): # check if we are at a terminal of a strand
                    consecutive_count += 1 # add terminal nucleotide
                    consecutive_stacked[scenario][t].append(consecutive_count)
                    prev  = -1
                    consecutive_count = 1
                    continue
                else:
                    consecutive_count += 1
                prev = coord_value
            
            n_broken_stacking[scenario].append(broken_count)

            if not consecutive_stacked[scenario][t]: # if there are NO stacked bases
                consecutive_stacked[scenario][t].append(0)

    return n_broken_stacking, consecutive_stacked

def plot_histogram(consecutive_stacked):
    fig      = plt.figure()
    ax       = fig.add_axes(111, projection='3d')
    binwidth = 1
    for i in range(len(consecutive_stacked)):
        data = sum(consecutive_stacked[i], [])
        n, bins = np.histogram(data,
                               density=True,
                               bins=np.arange(min(data)-binwidth/2, max(data)+1+binwidth/2, binwidth))
        bincenters = 0.5*(bins[1:] + bins[:-1])
        ax.bar3d(bincenters, i-0.2, 0, 0.3, 0.1, n)

    # set viewing angle
    ax.view_init(azim=-50) # Changes the elevation (elev) and azimuth (azim)

    # aspect ratio
    # x, y, z?
    ax.set_box_aspect((11, 8, 5))

    # make background transparent
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    # x-axis
    ax.set_xticks(np.arange(0, 22, 2), np.arange(0, 22, 2), rotation=0, ha='center', va='baseline')
    ax.set_xticks(np.arange(1, 22, 2), minor=True)
    ax.tick_params(axis='x', which='major', pad=0)
    ax.set_xlabel("Number of consecutively-stacked nucleotides", labelpad=-2)

    # y-axis
    ax.set_yticks(np.arange(len(legend)), legend, rotation=0, ha='left', va='baseline')
    ax.tick_params(axis='y', which='major', pad=-2)
    ax.set_ylabel("Scenarios", labelpad=-2)

    # z-axis
    ax.set_zticks(ax.get_zticks(), ax.get_zticklabels(), rotation=0, ha='center', va='baseline')
    ax.tick_params(axis='z', which='major', pad=-1)
    ax.set_zlabel("Probability density")
    ax.zaxis._axinfo['juggled'] = (1,2,0) # z-axis left side

    plt.savefig("stacking_hist.svg", bbox_inches="tight", dpi=600)

################################################################################################
#
# MAIN PROGRAM
#
################################################################################################

def main():
    # get data from .xvg files
    time, stacking_coords                  = get_data(paths, com_dir, vec_dir)
    n_broken_stacking, consecutive_stacked = analyze_data(stacking_coords)

    # set rcParams
    font_leg = set_rcParameters()
    
    # set figure dimensions
    fig_width   = 5
    golden_mean = (np.sqrt(5)-1.0)/2.0     # aesthetic ratio
    fig_height  = fig_width*golden_mean    # height in inches

    # plot the data
    plot_data(time,
              n_broken_stacking,
              "Simulation time (ns)",
              "Number of broken stacking",
              "Simulated Annealing",
              legend,
              "broken_stacking_vs_time.svg",
              fig_width,
              fig_height)

    plot_histogram(consecutive_stacked)
    
    # print statistics
    for i in range(len(n_broken_stacking)):
        print("Average number of broken stacking for file " + str(i+1) + ": " + str(round(statistics.mean(n_broken_stacking[i]),1)) + " +/- " + str(round(statistics.stdev(n_broken_stacking[i]),1)))
        
    for i in range(len(n_broken_stacking)):
        print("Average number of broken stacking for file " + str(i+1) + " (excluding first 600 ns): " + str(round(statistics.mean(n_broken_stacking[i][12000:]),1)) + " +/- " + str(round(statistics.stdev(n_broken_stacking[i][12000:]),1)))

if __name__ == "__main__": 
    main()
    