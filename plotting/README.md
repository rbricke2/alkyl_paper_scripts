# Overview 

A collection of scripts used to make plots for the paper.

# Programs

* `functions_for_plots.py`: Function file containing functions that multiple scripts use.
* `plot_experimental_melting_temp_data.py`: Plots Souyma Chandrasekhar's melting temperature data.
* `plot_hbond.py`: Plots 2D color plots showing the existence of Watson-Crick hydrogen bonding between base pairs throughout the duplex for each frame of the simulation. Needs `.xvg` files outputted by the GROMACS utilities `distance` and `angle`.
* `plot_radius_of_gyration.py`: Plots the radius of gyration as a function of time. Input file(s) assumed to be `.xvg` file(s) outputted by the GROMACS utility `gyrate`.
* `plot_stacking.py`:
* `plot_temp.py`:
