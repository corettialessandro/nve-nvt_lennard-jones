# nve-nvt_lennard-jones

Small code intended to perform molecular dynamics simulations of Lennard-Jones particles in NVE and NVT ensemble. 
Based on a source code by Prof. Mauro Ferrario (UNIMORE - Università degli Studi di Modena e Reggio Emilia)

It requires numba (https://numba.pydata.org) for the compilation of the force routine (calcener.py) in addition of the standard python routines for scientific computing (numpy, for example).

All routines needed in a MD code are contained in MDLJ_lib.py except for the one used to compute forces, which is contained in calcener.py.
The files NVE_lj.py and NVT_lj.py call the routines in MDLJ_lib.py and contain the routines for perfoming a full MD simulation in the respective ensemble.
Finally, the files NVE_main.py and NVT_main.py contain the call to the functions in the two previous files and set the parameters of the simulation.
