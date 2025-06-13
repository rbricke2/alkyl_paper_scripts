# Program

* `make_hbond_index_files.py`: Creates the index files `hbond_dist.ndx` and `hbond_angle.ndx` used as input for GROMACS utilities `distance` and `angle`, respectively. Only one donor-hydrogen-acceptor triplet, i.e. the one where the acceptor is a nitrogen atom, of a Watson-Crick base pair is considered. The file `hbond_dist.ndx` lists the atom IDs of the donor and acceptor atoms. The file `hbond_angle.ndx` lists the atom IDs of each triplet in this order: hydrogen, donor, then acceptor.
* `make_nucleobase_plane_COM_index_files.py` (not used for paper):  Creates index files `nucleobase_vec_atoms.ndx` and `nucleobase_COM_atoms.ndx` used as input for GROMACS utility 'traj'. The file `nucleobase_vec_atoms.ndx` lists the atoms IDs which are the endpoints of vectors $a$ and $b$ (refer to <cite>[this paper][1]</cite>) so that we can retrieve their positions in space. The file `nucleobase_COM_atoms.ndx` lists the ID of each heavy atom (i.e. non-hydrogen atom) in each nucleobase for center of mass calculation.

[1]: https://doi.org/10.1021/ct501025q

