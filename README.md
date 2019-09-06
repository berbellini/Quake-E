###################################
# READ ME for the quake-E package #
###################################

Created by Andrea Berbellini, September 2019
For questions, bugs, comments please email me: andrea.berbellini@ingv.it 

###  Introduction
The goal of this package is to measure Rayleigh wave ellipticity from earthquake data,
as in Berbellini et al. 2016.

If used, please cite us:
Berbellini, A., Morelli, A., Ferreira, A.M.G.; Ellipticity of Rayleigh waves in basin and hard-rock sites in Northern Italy, Geophysical Journal International, Volume 206, Issue 1, July, 2016, Pages 395–407, https://doi.org/10.1093/gji/ggw159


###  Prerequisites
- Python2.7
- Obspy1.1.0
- Matplotlib2.2.3

###  Installation
The code does not need to be installed or compiled, since it is completelly written in Python2.7

###  Run a measurement
The code comes with a dataset of 804 teleseisms recorded by station GIMEL (St. Georges, Gimel, Switzerland). 
Data have been provided by the Swiss Seismological Service via ORFEUS. 
The seismograms are ready-to-use, already corrected for the instrument response and rotated toward the event 
epicentre. Quake-E uses only the vertical and radial components.
For reference, precomputed results are stored in the folder “Results_DEMO”

To run a measurement from scratch please follow these steps:
1. Run the main routine typing:
> python start_measure.py

The code will start measuring ellipticity from all the available records, looping along all the chosen periods.

After approx 1 hour the raw results will be saved in the folder “Results/GIMEL/”, one file for each period.
Each resulting file contains the measurements for each earthquake using the format:
Event code | Ellipticity measurement | Standard deviation | Signal-to-noise ratio | Maximum value of the characteristic function.

2. After the computation is done, move to the “post_processing” folder. Here run the code to analyse the results and plot them.
> python process_and_plot.py

This code calculates the median ellipticity and related error for each period and saves the results 
in the folder “Results/GIMEL/” in the file named “ellipticity_curve.txt”

It also plots the results in 2 different ways:
A plot showing all the measurements done with errorbars (“all_measurements_TXXs.png”) and the whole 
ellipticity curve. It also compare the results obtained in this run with previous results from Berbellini et al. 2016. 

3. to plot the histograms of the measurements simply type:
> python plot_histograms.py
The resulting histogram will be stored in the same folder.

###  Changing dataset
In order to use your own measurements, store the seismograms in .sac format, rotated toward the event epicentre in order 
to obtain the vertical and radial component. Then edit the start_measure.py routine accordingly.



###  References
Berbellini, A., Morelli, A., Ferreira, A.M.G.; Ellipticity of Rayleigh waves in basin and hard-rock sites in Northern Italy, Geophysical Journal International, Volume 206, Issue 1, July, 2016, Pages 395–407, https://doi.org/10.1093/gji/ggw159

M. Beyreuther, R. Barsch, L. Krischer, T. Megies, Y. Behr and J. Wassermann (2010)
ObsPy: A Python Toolbox for Seismology SRL, 81(3), 530-533 DOI: 10.1785/gssrl.81.3.530

Swiss Seismological Service (SED) At ETH Zurich. (1983). National Seismic Networks of Switzerland. ETH Zürich. https://doi.org/10.12686/sed/networks/ch 
