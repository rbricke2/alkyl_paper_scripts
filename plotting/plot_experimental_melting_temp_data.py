"""
    Author: Rachel Bricker
    Year:   2025
"""

"""
    Plots Soumya's experimental data
"""

"""
   usage: python3 plot_experimental_melting_temp_data.py
"""

import numpy as np 
from functions_for_plots import plt, font_manager, set_rcParameters

################################################################################################
#
# MAIN PROGRAM
#
################################################################################################

def main():
    # rgb colors used for plots
    c_butyl  = (0.129, 0.588, 0.953)
    c_heptyl = (0.957, 0.263, 0.212)
    c_decyl  = (0.584, 0.647, 0.651)
    c_unmod  = (1,     0.498, 0.055)
    c_pentyl = (0.941, 0.894, 0.259)

    font_leg = set_rcParameters()
    
    fig_height = 1.8  
    fig_width  = 3.2
    
    ################################################
    # plot bar graph
    
    # data 
    n_groups = np.arange(3) 
    group_1  = [55, 55, 59] 
    group_2  = [51, 54, 57] 
    group_3  = [50, 51, 55] 
    width    = 0.2
    
    # initialize figure
    fig, ax = plt.subplots(1, figsize=(fig_width, fig_height))
    
    # plot data in grouped manner of bar type
    plt.bar(n_groups-width, group_1, width, color=c_butyl,  edgecolor="black", linewidth=0.8)
    plt.bar(n_groups,       group_2, width, color=c_heptyl, edgecolor="black", linewidth=0.8)
    plt.bar(n_groups+width, group_3, width, color=c_decyl,  edgecolor="black", linewidth=0.8)
    plt.xticks(n_groups, ['4PS\nspaced', '2PS\ntogether', '4PS\ntogether']) 
    plt.xlabel("Position of phosphorothioate") 
    plt.ylabel("Melting Temp. (${}^{\circ}$C)") 

    plt.legend(["Butyl", "Heptyl", "Decyl"], loc='center left', bbox_to_anchor=(1, 0.5), prop=font_leg)

    plt.grid(axis='y')

    plt.tight_layout()

    # save figure
    plt.savefig("soumya_barplot.svg", bbox_inches="tight", dpi = 600)
    
    ################################################
    # plot scatter plot
    
    # data
    x  = [   0,  1,  2,  4,  6,  8, 10]
    y1 = [  70, 69, 68, 66, 65, 63, 62]
    y2 = [None, 61, 58, 55, 52, 47, 43]
    y3 = [None, 61, 58, 54, 51, 45, 42]
    y4 = [None, 61, 58, 53, 45, 43, 40]

    # initialize figure
    fig, ax = plt.subplots(1, figsize=(fig_width, fig_height))

    # plot data
    plt.scatter(x, y1, marker=",", facecolor=c_unmod,  edgecolors='black', linewidth=0.8, s=60)
    plt.scatter(x, y2, marker="o", facecolor=c_butyl,  edgecolors='black', linewidth=0.8, s=60)
    plt.scatter(x, y3, marker="X", facecolor=c_pentyl, edgecolors='black', linewidth=0.8, s=60)
    plt.scatter(x, y4, marker="^", facecolor=c_heptyl, edgecolors='black', linewidth=0.8, s=60)
    plt.xlabel("Number of phosphorothioates") 
    plt.ylabel("Melting Temp. (${}^{\circ}$C)") 

    plt.legend(["Unmodified", "Butyl", "Pentyl", "Heptyl"], loc='center left', bbox_to_anchor=(1, 0.5), prop=font_leg) 

    ax.minorticks_on()

    plt.grid()

    plt.tight_layout()

    # save figure
    plt.savefig("soumya_scatterplot.svg", bbox_inches="tight", dpi = 600)
    

if __name__ == "__main__": 
    main()

