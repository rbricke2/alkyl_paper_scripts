# Function file for plotting graphs
# usage: from functions_for_plots import *
# Author: Rachel Bricker

import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

def set_rcParameters():
    font_size    = 8
    font_leg     = font_manager.FontProperties(family="Arial", style='normal', size=font_size)
    widths       = 1
    line_width   = widths+(widths*0.5)
    rcParameters = {'font.size'        : font_size,     # font size
                    'font.family'      : "sans-serif",  # specify the font
                    'font.sans-serif'  : "Arial",
                    'axes.linewidth'   : widths,        # set width of axes
                    'xtick.major.width': widths,        # set tick widths
                    'xtick.minor.width': widths,
                    'ytick.major.width': widths,
                    'ytick.minor.width': widths,
                    'grid.color'       : 'lightgray',   # set grid color
                    'grid.linewidth'   : widths,        # set grid width
                    'grid.linestyle'   : '--',
                    'lines.linewidth'  : line_width,    # line width
                    'patch.linewidth'  : widths,        # legend line width
                    'axes.axisbelow'   : True}          # plot grid behind data
    plt.rcParams.update(rcParameters)
    
    return font_leg
