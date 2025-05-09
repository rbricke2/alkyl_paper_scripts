"""
    Author:          Rachel Bricker
    Year:            2025
    Acknowledgments: Thanks to Dr. John Portman for helping me.
"""

r"""
    Plots stacking analysis.

    Base stacking is characterized geometrically. The stacking definition proposed by the Florian group [1] is used.
    That is, the distance between the centers of mass of two nucleobases and the angle between their planes are combined to produce
    a single reaction coordinate, $\xi$. Adjacent nucleobases are considered stacked if $\xi$ is less than or equal to 0.6 nm.

    The plane of a nucleotide was determined using the definition presented by the Turner group [2]. That is, the plane of a nucleotide
    is defined by two vectors, $a$ and $b$, whose cross products $a \times b$ and $b \times a$ define each base’s normal vectors. For purines,
    $a$ is defined from the center of mass to the C8 atom and $b$ is defined from the center of mass to the nitrogen (resp. oxygen) attached to
    the C6 atom in adenine (resp. guanine). For pyrimidines, $a$ is defined from the center of mass to the oxygen group attached to the C2 atom
    and $b$ is defined from the center of mass to the oxygen (resp. nitrogen) attached to the C4 atom in thymine (resp. cytosine). The atoms are
    surrounded by asterisk symbols in the below images.

            THYMINE                                          CYTOSINE
                                                                        H42  H41
                                                                          \  /
                  H51    *O4*                                             *N4*           
                   |      ||                                               | 
               H52-C5M    C4     H3                                        C4  
                   |  \ /    \  /                                        /    \\    
                  H53  C5     N3                                     H5-C5     N3   
                       ||     |                                         ||     |  
                    H6-C6     C2                                     H6-C6     C2
                        \    / \\                                        \    / \\
                          N1   *O2*                                        N1   *O2* 
                           \                                                \        
                            ~                                                ~



            ADENINE                                            GUANINE
                   H61  H62                                          *O6*                  
                     \  /                                             ||                     
                     *N6*                                             C6                     
                      |                                             /    \                  
                      C6                                        H1-N1     C5--N7 \\       
                   //    \                                         |      ||     *C8*-H8
                  N1      C5--N7 \\                                C2     C4--N9 /         
                  |       ||     *C8*-H8                          / \\   /      \           
                  C2      C4--N9 /                          H21-N2     N3        ~          
                 / \\    /     \                                |                           
               H2     N3        ~                              H22                            
                                 

    The center of mass of each nucleobase for each snapshot was obtained using the standard GROMACS utility ‘traj’; hydrogen atoms were not included
    in the calculation of the center of mass. The angle between planes of two nucleobases and the distance between their mass-centers is measured
    using this script written in-house.

    [1] Jafilan, S., Klein, L., Hyun, C., & Florian, J. (2012). Intramolecular base stacking of dinucleoside monophosphate anions in aqueous solution.
        The Journal of Physical Chemistry B, 116(11), 3613-3618.
    [2] Condon, D. E., Kennedy, S. D., Mort, B. C., Kierzek, R., Yildirim, I., & Turner, D. H. (2015). Stacking in RNA: NMR of four tetramers benchmark molecular dynamics.
        Journal of chemical theory and computation, 11(6), 2729-2742.
"""

"""
   usage: python3 plot_stacking.py
      1. list of model names for legend (must be parallel w.r.t the arguments i to i+n)
      2. name of directory holding .xvg files containing the COM of each nucleobase outputted
         by GROMACS utility `traj`
      3. name of directory holding .xvg files containing the x, y, z position of the atoms that define
         the vectors in each nucleobase outputted by GROMACS utility `traj`
      4. total number of nucleotides
      5. enter 1 if double-stranded, 0 if single-stranded
      i. path to directories that contain arguments (2.) and (3.) that you want to plot

   example: python3 plot_stacking.py \
            "(a),(b),(d)" \
            com_files \
            vec_files \
            42 \
            1 \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/annealing_AMBER/dsDNA1/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/annealing_AMBER/dsDNA3/ \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/annealing_AMBER/dsDNA4/ 
"""

import sys
import numpy as np 
from functions_for_plots import *

# command line input
input_list = sys.argv[1].replace("\\n", "\n").replace("\\t", "\t").replace("\\(", "(").replace("\\)", ")")
legend     = input_list.split(',')
com_dir    = str(sys.argv[2])
vec_dir    = str(sys.argv[3])
n_residues = int(sys.argv[4])
ds         = bool(int(sys.argv[5]))
paths      = list(sys.argv[6:])

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

def get_data(paths, n_residues, com_dir, vec_dir, ds):
    stacking_coords = [ [] for scenario in range(len(paths)) ]
    for scenario in range(len(paths)):
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
    ax.set_ylabel("Models", labelpad=-2)

    # z-axis
    ax.set_zticks(ax.get_zticks(), ax.get_zticklabels(), rotation=0, ha='center', va='baseline')
    ax.tick_params(axis='z', which='major', pad=-1)
    ax.zaxis._axinfo['juggled'] = (1,2,0) # z-axis left side
    ax.set_zlabel("Probability density", labelpad=-2)

    plt.savefig("stacking_hist.svg", bbox_inches="tight", dpi=600)

################################################################################################
#
# MAIN PROGRAM
#
################################################################################################

def main():
    # get data from .xvg files
    time, stacking_coords                  = get_data(paths, n_residues, com_dir, vec_dir, ds)
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
              "Number of broken stacking interactions",
              "Simulated Annealing",
              legend,
              "broken_stacking_vs_time.svg",
              fig_width,
              fig_height)

    plot_histogram(consecutive_stacked)

    # print statistics
    frame_600ns = 12000
    for i in range(len(n_broken_stacking)):
        print("Average number of broken stacking for file " + str(i+1) + ": " + str(round(statistics.mean(n_broken_stacking[i]),1)) + " +/- " + str(round(statistics.stdev(n_broken_stacking[i]),1)))
    for i in range(len(n_broken_stacking)):
        if len(n_broken_stacking[i]) > frame_600ns+1:
            print("Average number of broken stacking for file " + str(i+1) + " (excluding first 600 ns): " + str(round(statistics.mean(n_broken_stacking[i][frame_600ns:]),1)) + " +/- " + str(round(statistics.stdev(n_broken_stacking[i][frame_600ns:]),1)))

if __name__ == "__main__": 
    main()
    