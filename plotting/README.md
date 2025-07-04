# Overview 

A collection of scripts used to make plots for the paper.

# Programs

* `functions_for_plots.py`: Function file containing functions that multiple scripts use.
* `plot_experimental_melting_temp_data.py`: Plots Souyma Chandrasekhar's melting temperature data.
* `plot_hbond.py`: Plots 2D color plots showing the existence of Watson-Crick hydrogen bonding between base pairs throughout the duplex for each frame of the simulation. Needs `.xvg` files outputted by the GROMACS utilities `distance` and `angle`.
* `plot_radius_of_gyration.py`: Plots the radius of gyration as a function of time. Input file(s) assumed to be `.xvg` file(s) outputted by the GROMACS utility `gyrate`.
* `plot_stacking.py` (not used in paper): Analyzes stacking between each base pair step. The stacking definition proposed by the <cite>[Florian group][3]</cite> is used. The plane of a nucleotide was determined using the definition presented by the <cite>[Turner group][4]</cite>.
* `plot_x3DNA.py`: Plots twist averaged over the base pairs of DNA, excluding the three terminal ones at each end of the duplex. The twist of each base pair step is calculated using the <cite>[3DNA][1]</cite> and <cite>[do_x3dna][2]</cite> softwares.

[1]: https://doi.org/10.1093/nar/gkg680
[2]: https://doi.org/10.1093/bioinformatics/btv190
[3]: https://doi.org/10.1021/jp209986y
[4]: https://doi.org/10.1021/ct501025q
