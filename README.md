# plosCB-ponzi-HD
Code and data for: Ponzi et al. PlosCB, Striatal Network Modeling in Huntington’s Disease (2020)

Data files.

Spike data generated from the model used in the paper.

Files are named, spike_output_gEx_gIy.mat, where x and y denote the values of the corresponding parameters gE and gI used in the paper. Due to their large size a few files have been divided into parts, called spike_output_gEx_gIy_pn.mat, n=1,2,.. These are files: spike_output_gE110_gI1.mat, spike_output_gE75_gI1.mat, spike_output_gE50_gI1.mat and spike_output_gE40_gI1.mat.

There are two columns in each of these Matlab data files. The first column is the spike time, in units of 0.2 ms. The second column is the cell index number in the network (0 to 2499). Simulation time length is 212 seconds. Results in the paper discarded the first 10 seconds.

Code files.

Code used to generate spike data used in the paper in MGS.

Mouse data files.

Experimental mouse recording data used in the paper in Matlab files.  
