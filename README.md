[![Paper](https://img.shields.io/badge/paper-arXiv%3A2105.06402-B31B1B.svg)](https://arxiv.org/abs/2105.06402)
[![DOI](https://zenodo.org/badge/507527663.svg)](https://zenodo.org/badge/latestdoi/507527663)

# Essential role of magnetic frustration in the phase diagrams of doped cobaltites

[Peter P. Orth](https://faculty.sites.iastate.edu/porth/), D. Phelan, J. Zhao, H. Zheng, J. F. Mitchell, C. Leighton, Rafael M. Fernandes

### Abstract 
Doped perovskite cobaltites (e.g., La$_{1-x}$ Sr$_x$ CoO$_3$) have been extensively studied for their spin-state physics, electronic inhomogeneity, and insulator-metal transitions. Ferromagnetically-interacting spin-state polarons emerge at low $x$ in the phase diagram of these compounds, eventually yielding long-range ferromagnetism. The onset of long-range ferromagnetism ($x \approx 0.18$) is substantially delayed relative to polaron percolation ($x \approx 0.05$), however, generating a troubling inconsistency. Here, Monte-Carlo simulations of a disordered classical spin model are used to establish that previously ignored _magnetic frustration_ is responsible for this effect, enabling faithful reproduction of the magnetic phase diagram.

### Description
This repository includes information, code, and data to generate the theoretical figures of the paper (Figs.3 and 4).

### Requirements
* g++
* openmp
* matplotlib

### Data Generation
* ```mcdis-pt.cpp``` C++ file to perform classical Monte Carlo (MC) simulation of disordered 3D Ising model on cubic lattice. MC updates are using a combination of local Metropolis updates and a parallel-tempering update. The program runs ```T_steps``` simulations in parallel in a specified temperature interval. Each simulation is repeated ```disorder_steps``` number of times with random quenched disorder. 

After compilation the file can be run with the command
```
./mcdis-pt dim L Lz J33 J34 J44 T_min T_max T_steps disorder_filling disorder_steps mc_bins mc_binsize thermalization_bins seed_disorder seed_mc
```
The parameters denote:
- ```dim```: dimensionality of the system (```dim=3``` for the cubic lattice). 
- ```L```: linear number of sites along the $x$ and the $y$ direction
- ```Lz```: linear number of sites along the $z$ direction
- ```J33, J34, J44```: spin exchange coupling constants
- ```T_min, T_max```: specifies range of temperatures in which system is simulated
- ```T_steps```: number of temperatures that are being simulated. These are spaced geometrically between ```T_min``` and ```T_max```
- ```disorder_filling```: fraction of core Co$^{4+}$ sites (site label $4$) that are populated randomly for a given quenched disorder realization
- ```disorder_steps```: number of independent quenched disorder realizations that are simulated. Output observables and bond/site fractions are averaged over these independent disorder realizations
- ```mc_bins```: number of MC bins during which observables are measured
- ```mc_binsize```: number of MC steps per MC bin. Each step consists of a single Metropolis MC step ($N$ update attempts) and a single parallel tempering step
- ```thermalization_bins```: number of mc_bins used for thermalization at the beginning. Observables are not measured during thermalization. The total number of MC bins is therefore ```mc_bins + thermalization_bins```
- ```seed_disorder```: integer seed of random number generator that sets quenched disorder configuration
- ```seed_mc```: integer seed of random number generator for MC simulations

### Data description
The program ```mcdis-pt``` outputs the following data files
* ```Data-AveragedObservables-L_#-Lz_#-disorder_filling_#-J33_#-J34_#-J44_#-TMin_#-TMax_#-TSteps_#.dat```: contains MC observables that averaged over quenched disorder configurations. See beginning of the file for a detailed description of the content of the different columns. Standard deviations are computed using the jackknife method. 

The program generates also intermediate output after every quenched disorder realization. There are two files that are output that are identical at the end of the calculation, but are different at temporary outputs to prevent problems in case the program is cancelled during writing output. 
* ```Data-Fraction.dat```: contains disorder averaged fraction of filled sites and bonds. See .py plotting scripts for more details.

Examples of these two files together with plotting scripts to plot select quantities from these files are located in subfolder ```figures```.

### Data plotting
Python matplotlib scripts for plotting the output data. These are used for generating the figures in the paper. These files together with the necessary data files are located in subfolder ```figures```.
* ```plot_fraction.py```: plots fraction of filled spinful sites for a given fraction of $4$ sites $x$. The inset uses data from file ```Data-Magnetization-L_10-Lz_10-J33_0.000000-J34_-0.500000-J44_-1.000000-temperature_0.100000-run_steps_1000.dat``` that contains the number of isolated polarons in $L = L_z = 10$ sized system, averaged over $1000$ quenched disorder realizatoins.
* ```plot_fraction_bonds.py```: plots fraction of filled bonds for a given fraction of $4$ sites $x$
* ```plot_phase_diagram.py```: plots phase diagram (using data file ```Data-Tc-J33_0.2-Mod.dat``` that is extracted from crossing of ferromagnetic Binder cumulant obtained in MC simulation). Plots also spin observables magnetization and susceptibilities.
