#############################################
#         				    				#
#   Andrea Berbellini - PhD student         #
#       INGV-Bologna   September 2019 	    #
#	andrea.berbellini@ingv.it           	#
#					    					#
#############################################
from E_function import *
import numpy as np
import glob
import os
import sys



channel = "HH"
station = "GIMEL"
network = "CH"

k = 1.4

periods = [10, 12.5, 15, 20, 25, 30, 35, 40, 50, 60, 80 ,100]


cmtpath = "./CMTSOLUTIONS"
results_folder = "Results/"
database_path = "./Data/" + station + "/"

os.system("touch " + results_folder + "log_file.txt")

log = open(results_folder + "log_file.txt", "a")
print "Measuring " + station
# creating output folder
os.system("mkdir "+ results_folder + station)
folder_list = os.listdir(database_path)
for T0 in periods:
	out = open(results_folder +  station + "/results_T"+str(T0)+".txt", "w")
	out.write("event" + "\t" + "E" + "\t" + "sd" + "\t" + "signal/noise" +"\t"+ "char_func_max" + "\n")
	n = 0
	m = 0
	for i in range(0,5):#len(folder_list)):
		n += 1
		event = folder_list[i]
		print "Event: ", event, i+1, "/", len(folder_list)
		Z_path = glob.glob(database_path + folder_list[i] + "/" + network + ".*Z*.SAC.corr")
		R_path = glob.glob(database_path + folder_list[i] + "/" + network + ".*R*.SAC.corr")

		if Z_path != [] and R_path != []:
			#try:
				T, E, sd, snr, char_func_max = E_function(R_path[0], Z_path[0], T0, "no", event, station, cmtpath)
			#	sys.exit()
				out.write(event + "\t" + str(E) + "\t" + str(sd) + "\t" + str(snr)+ "\t" + str(char_func_max) + "\n")
				print "Period: ", T0, "s Ellipticity: ", E
				print "----------------------"
				if  isinstance(E, float):
					m += 1
			#except:
			#	E = "Error"
			#	sd = "Error"
			#	snr = "Error"
			#	char_func_max = "Error"
			#	out.write(event + "\t" + E + "\t" + sd + "\t" + snr +"\t"+ char_func_max + "\n")
			#	pass
	
print "\n\n\n\n--------------------------------"
print "Tot. events number: ", n
print "Measurements performed: ", m
perc = 100 * m / n
print "Percentage: ", perc, "%"
print "--------------------------------\n\n\n\n"
log.write("Station: " + station + " Done correctly!\n")
#except:
#		log.write("Station: " + station + " had some problem :( \n")
#		pass


















	




