# RS-FS-LTS-clustering

This repository contains the code used for the manuscript:
Ter Wal, M., Tiesinga, P.H.E.,
Comprehensive characterization of oscillatory signatures in a model circuit with PV- and SOM-expressing interneurons (under review).

## Prerequisites
MATLAB (The MathWorks) with the Signal Processing Toolbox - The code has been developed on Matlab version R2018a.

For the Violin plots in Figure 5, the following function was used: 
violin.m 
Hoffmann H, 2015: violin.m - Simple violin plot using matlab default kernel density estimation. INRES (University of Bonn), Katzenburgweg 5, 53115 Germany.
https://uk.mathworks.com/matlabcentral/fileexchange/45134-violin-plot

## Getting started
Download the .zip and unzip, or clone to your favourite path.
Make sure you Matlab path is set to include the code.

The following two scripts do the bulk of the work:
- RunSimulations: runs the network simulation or simulations as defined in the script GetSettings. It also calls all the relevant functions to analyse the spiking and LFP behavior of the circuit. For each of the 20 circuit motifs used in the manuscript, a separate version of GetSettings is provided in the repo, to make replication easier.
- ClusterAnalysis: Once all simulations have been performed and saved, this script walks through the steps to perform k-means clustering and plots the results.

## License and disclaimer
This code is released under a CC BY 4.0 license, meaning that this code is free to use and modify for everyone, given appropriate credit, but comes without warranty.


Correspondence: Marije ter Wal - m.j.terwal@bham.ac.uk
