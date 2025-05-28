"""
    Author: Rachel Bricker
    Year:   2025
"""

"""
    Plots experimental data obtained by Soumya Chandrasekhar.
"""

"""
   usage: python3 plot_experimental_melting_temp_data.py
"""

import numpy as np
import csv
from functions_for_plots import plt, font_manager, set_rcParameters, statistics

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

    # set rcParams
    font_leg = set_rcParameters()
    
    # set figure dimensions
    fig_height = 1.8  
    fig_width  = 3.8
    
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
    
    # label axes
    plt.xlabel("Modification pattern") 
    plt.ylabel("$T_m$ (${}^{\circ}$C)") 

    # give legend information
    legend = plt.legend(["Butyl", "Heptyl", "Decyl"], loc='center left', bbox_to_anchor=(1, 0.5), prop=font_leg, frameon=False)

    # only show grid for y-axis
    plt.grid(axis='y')
    
    # change y-axis range
    plt.ylim((30,60))

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
    
    # label axes
    plt.xlabel("Number of phosphorothioates") 
    plt.ylabel("$T_m$ (${}^{\circ}$C)") 

    # give legend information
    plt.legend(["Unmodified", "Butyl", "Pentyl", "Heptyl"], loc='center left', bbox_to_anchor=(1, 0.5), prop=font_leg, frameon=False) 

    # show minor ticks
    ax.minorticks_on()

    # show grid
    plt.grid()

    plt.tight_layout()

    # save figure
    plt.savefig("soumya_scatterplot.svg", bbox_inches="tight", dpi = 600)

    ################################################
    # plot derivative curves

    with open('PTO_triplicates.csv', mode='r') as file:
        csv_file = csv.reader(file)
        rows     = list(csv_file)

    # remove first and last rows
    rows = rows[1:-1]

    # retrieve data
    x           = [int(item[0]) for item in rows]
    y_no_PS_avg = [statistics.mean( [float(item[14]), float(item[15]), float(item[16])]) for item in rows]
    y_no_PS_std = [statistics.stdev([float(item[14]), float(item[15]), float(item[16])]) for item in rows]
    y_1_PS_avg  = [statistics.mean( [float(item[17]), float(item[18]), float(item[19])]) for item in rows]
    y_1_PS_std  = [statistics.stdev([float(item[17]), float(item[18]), float(item[19])]) for item in rows]
    y_2_PS_avg  = [statistics.mean( [float(item[20]), float(item[21]), float(item[22])]) for item in rows]
    y_2_PS_std  = [statistics.stdev([float(item[20]), float(item[21]), float(item[22])]) for item in rows]
    y_4_PS_avg  = [statistics.mean( [float(item[23]), float(item[24]), float(item[25])]) for item in rows]
    y_4_PS_std  = [statistics.stdev([float(item[23]), float(item[24]), float(item[25])]) for item in rows]

    y = [y_no_PS_avg, y_1_PS_avg, y_2_PS_avg, y_4_PS_avg]
    e = [y_no_PS_std, y_1_PS_std, y_2_PS_std, y_4_PS_std]
    
    #print(y[2][16])
    #print(e[2][16])

    del y_no_PS_avg, y_1_PS_avg, y_2_PS_avg, y_4_PS_avg, y_no_PS_std, y_1_PS_std, y_2_PS_std, y_4_PS_std

    # figure dimensions
    fig_width  = 7.2
    fig_height = 3.6

    # set rcParams
    font_leg = set_rcParameters()

    # set figure dimensions
    fig, axes = plt.subplots(1, 4, sharey=True, figsize=(fig_width, fig_height))

    # titles
    titles = ["Unmodified", "One\nPhosphorothioate", "Two\nPhosphorothioates", "Four\nPhosphorothioates"]

    # plot data
    for i in range(len(axes)):
        # line
        axes[i].plot(x, y[i], color="r", alpha=0.3, linewidth=1)
        # error
        axes[i].errorbar(x, y[i], yerr=e[i], ls = "None", capsize=3, color='black', elinewidth=0.5, markeredgewidth=0.5)
        # scatter
        axes[i].scatter(x, y[i], s=4, facecolor="r", zorder=2, edgecolors="black", linewidth=0.5)

        # set y-axis label on leftmost plot
        if i == 0:
            axes[i].set_ylabel("Average Derivative")

        # show grid
        axes[i].grid()

        # set title
        axes[i].set_title(titles[i])

    # add a big axis and hide frame for common x-axis label
    fig.add_subplot(111, frameon=False)

    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)

    # common x-axis label
    plt.xlabel("Cycles")

    plt.tight_layout()

    # save figure
    plt.savefig("derivative.svg", bbox_inches="tight", dpi=600)

if __name__ == "__main__": 
    main()

