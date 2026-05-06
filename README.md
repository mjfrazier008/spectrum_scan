# spectrum scan
Spectral calculations for continuous photonic systems as described in [Effects of interface regularity on the bulk-edge correspondence in continuum photonic systems](https://arxiv.org/abs/2605.02100). 

## spectrum_calc_scan.m
Use this script to calculate spectra. Requires MATLAB parallel computing toolbox as written. To remove this toolbox requirement change the annotated lines in "residue_map.m" and "grid_adapt.m". 

## cold_plasma_tm_(...).m
Functions which define the appropriate interface operator (operator $V(k_x, E)$ as described in the above reference). Suffixes correspond to various models used. Change parameter values $\omega_p, k_c, \beta$ inside these functions. 

## H_tm_hydro.m
Use to generate and diagonalize bulk (constant coefficient) Hamiltonians of the hydrodynamic system.

## test_modes_vec.m
Can be used to extract eigenvectors at a particular $(k_x, E)$ value.

## scatter_sort.m
Used to sort output of spectral calculations into data which can be plotted as a scatter plot.
