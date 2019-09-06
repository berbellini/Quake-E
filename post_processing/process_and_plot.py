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

outfile = open("../Results/"+station+"/ellipticity_curve.txt","w" )
outfile.write("Period\tE\tError\n")
EEmedian, ee_error = [],[]
for T in periods:
	filename = results_folder + "results_T"+str(T)+".txt"
	lines = open(filename).readlines()
	n = 0
	EE, ssd, nn = [],[],[]
	for i in range(1,len(lines)):
		try:
			E = float(lines[i].split()[1])
			E = np.log10(E)
			sd = float(lines[i].split()[2])
			snr = float(lines[i].split()[3])

			if snr > 100 and abs(E) < 1.:
				EE.append(E)
				ssd.append(sd)
				nn.append(n)
			n+=1
		except ValueError:
			pass

	Emedian = round(np.median(EE),5)
	Epercentile1 = np.percentile(EE, 15.9)
	Epercentile2 = np.percentile(EE, 84.1)
	err_inf = abs(Epercentile1 - Emedian)
	err_sup = abs(Epercentile1 - Emedian)
	error = round((err_inf + err_sup) / 2.0,5)


	EEmedian.append(Emedian)
	ee_error.append(error)

	#=========== write ellipticity curve to file =================
	outfile.write(str(T)+"\t"+str(Emedian)+"\t"+str(error)+"\n")

	plt.scatter(nn, EE, color="C0",s=5)
	plt.errorbar(nn, EE, yerr=ssd, fmt=" ", color="C0",alpha=0.2)
	plt.axhline(Emedian, color="red", label="median")
	plt.axhline(Epercentile1, color="blue", label="percentile 15.9")
	plt.axhline(Epercentile2, color="blue", label= "percentile 84.1")
	plt.suptitle("Station: "+station+ " Period: " + str(T)+ "s")
	plt.legend()
	plt.ylabel("Log(E)")
	plt.xlabel("# Event")

	plt.ylim(-1,1)

	plt.savefig("../Results/all_measurements_T"+str(T)+"s.png")
	plt.close()



outfile.close()



#=========== plotting ellipticity curve ###################
fig, ax = plt.subplots()
plt.scatter(periods, EEmedian, color="black", label="E (new measurements)")
plt.errorbar(periods, EEmedian, yerr=ee_error, color="black", fmt=" ")



lines = open("../GIMEL_Berbellini_et_al_2016.txt").readlines()
TT, EE = [],[]
for i in range(1,len(lines)):
	if float(lines[i].split()[0]) < 120:
		TT.append(float(lines[i].split()[0]))
		EE.append(float(lines[i].split()[1]))

plt.plot(TT, EE,color="C1", label="E from Berbellini et al. 2016",zorder=0)
plt.legend()
plt.ylabel("Log(E)")
plt.xlabel("Period (s)")

plt.ylim(-0.5, 0.5)
plt.xscale("log")
plt.suptitle("Station: "+station)
ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax.xaxis.set_minor_formatter(StrMethodFormatter('{x:.0f}'))
plt.savefig("../Results/Ellipticity_curve.png")
plt.close()






