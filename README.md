# plosCB-ponzi-HD
Code and data for: Ponzi et al. PlosCB, Striatal Network Modeling in Huntington’s Disease (2020)

Data files.

Spike data generated from the model used in the paper.

Spike_Output.dat, in each of six included directories, gE_15, gE_25, gE_40, gE_50, gE_75, and gE_110, include spike times from cell models in network models simulated at the corresponding parameter value, g_E.

There are two columns in each of these text data files. The first column is the is cell indices in the network (0 to 2499). The second column is the spike times for each cell, in units of 0.2 ms.  

In each file are multiple simulation outputs, concatenated, at increasing levels of inhibition (g_I). The first simulation is at g_I=1, and subsequent simulations then increase in steps of 4, i.e. g_I={1, 5, 9, ... }. 

New simulations can be determined by the fact that the spike time (first column) decreases. Simulation time length is 400 seconds.

Code files.

Code used to generate spike data used in the paper in MGS.

Mouse data files.

Experimental mouse recording data used in the paper in Matlab files.  
