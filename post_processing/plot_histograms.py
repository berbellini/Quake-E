#############################################
#         				    				#
#   Andrea Berbellini - PhD student         #
#       INGV-Bologna   September 2019 	    #
#	andrea.berbellini@ingv.it           	#
#					    					#
#############################################
import os
import sys
import numpy as np
import matplotlib.pylab as plt
import glob
from matplotlib.ticker import StrMethodFormatter, NullFormatter



station = "GIMEL"
results_folder = "../Results/" + station + "/"

periods = [10, 12.5, 15, 20, 25, 30, 35, 40, 50, 60, 80 ,100]
n = 1


# Read medians 
lines = open("../Results/"+station+"/ellipticity_curve.txt").readlines()
EEmedian = []
for i in range(1,len(lines)):
	Emedian = float(lines[i].split()[1])
	EEmedian.append(Emedian)


plt.figure(1,figsize=(8.27,11.69))
plt.subplots_adjust(left=None, bottom=None, right=None, top=0.95, wspace=None, hspace=0.5)
for T in periods:
	filename = results_folder + "results_T"+str(T)+".txt"
	lines = open(filename).readlines()
	
	EE, ssd, nn = [],[],[]
	for i in range(1,len(lines)):
		try:
			E = float(lines[i].split()[1])
			E = np.log10(E)
			sd = float(lines[i].split()[2])
			snr = float(lines[i].split()[3])
			if snr > 100 and abs(E) < 1.:
				EE.append(E)
		except ValueError:
			pass

	plt.subplot(len(periods),1,n)
	bins = np.arange(-1,1,0.02)
	values = plt.hist(EE, bins=bins)
	plt.axvline(EEmedian[n-1],color="C6",linewidth=3)
	max_freq = max(values[0])

	plt.text(-1, .5*max_freq, str(T)+"s", fontsize=15)



	n+=1

plt.subplot(len(periods),1,len(periods))
plt.xlabel("log(E)",fontsize=15)
plt.savefig("../Results/histograms.png")
plt.close()