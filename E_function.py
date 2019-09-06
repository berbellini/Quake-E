#!/usr/bin/python
#############################################
#         				    				#
#   Andrea Berbellini - PhD student         #
#       INGV-Bologna   September 2019 	    #
#	andrea.berbellini@ingv.it           	#
#					    					#
#############################################
import matplotlib.pyplot as plt
import numpy as np
from obspy.core import read
import obspy.signal
from scipy.interpolate import interp1d
import sys
from functions import *
from obspy.signal.cross_correlation import xcorr

def E_function(R_file, Z_file, T0, plot, event, station, cmtpath):
	threshold = 0.6
	k=1.4
	seconds_of_noise = 2400
	Tmin = T0 / np.sqrt(k)
        Tmax = T0 * np.sqrt(k)
	Fmin = 1/Tmax
	Fmax = 1/Tmin

	ell_sum = 0.0
	num_data = 0.0
	av_ell = 0.0
	sd_total = 0.0

	R_st = read(R_file)
	Z_st = read(Z_file)
	R_tr = R_st[0]
	Z_tr = Z_st[0]
	channel = Z_tr.stats.channel
	#-----------------------------
	# Preliminary check: delta
	#-----------------------------
	if R_tr.stats.delta !=  Z_tr.stats.delta:
		#print "Traces have different delta!"
		E = "DiffDelta"
		sd = "DiffDelta"
		snr = "DiffDelta"
		char_func_max = "DiffDelta"
		return T0, E, sd, snr,  char_func_max
	
	delta = Z_tr.stats.delta
	#-------------------------------------------------------
	# trimming both the trace to have the same npts	
	#-------------------------------------------------------
	start_time = Z_tr.stats.starttime
	end_time = Z_tr.stats.endtime
	R_tr.trim(start_time+100, end_time-100)
	Z_tr.trim(start_time+100, end_time-100)

	#----------------------------------
	# check npts
	#----------------------------------
	if R_tr.stats.npts != Z_tr.stats.npts:
		#print "Traces have different npts!"
		E = "DiffNPTS"
		sd = "DiffNPTS"
		snr = "DiffDelta"
		char_func_max = "DiffNPTS"
		return T0, E, sd, snr, char_func_max
	npts = Z_tr.stats.npts

	#--------------------------------------------------------
	# filtering traces
	#--------------------------------------------------------
	
	R_tr.detrend('demean')
	Z_tr.detrend('demean')
	R_tr.taper(0.05,type='cosine')
	Z_tr.taper(0.05,type='cosine')
	R_tr.filter('bandpass', freqmin = Fmin, freqmax = Fmax, corners=2, zerophase=True)
	Z_tr.filter('bandpass', freqmin = Fmin, freqmax = Fmax, corners=2, zerophase=True)
				
	#-----------------------------------------------
	# Signal/noise check
	#-----------------------------------------------
	R_noise = sig_noise(R_tr, seconds_of_noise)
	Z_noise = sig_noise(Z_tr, seconds_of_noise)
	snr = (R_noise + Z_noise) / 2.0

	#--------------------------------------
	# Start measurement
	#--------------------------------------				
	Z_dephase = Z_tr.copy()
	Z_dephase = phase_shift(Z_tr, T0/4.0)
	
	T_window = 3*T0
	cross_corr_trace = cross_correlation(R_tr, Z_dephase, T_window)

	R_env = R_tr.copy()
	R_env.data = obspy.signal.filter.envelope(R_env.data)
	Z_env = Z_dephase.copy()
	Z_env.data = obspy.signal.filter.envelope(Z_dephase.data)
	R_env_norm = R_env.copy()
	R_env_norm.normalize()
	Z_env_norm = Z_env.copy()
	Z_env_norm.normalize()

	ZR_env_norm = Z_env_norm.copy()
	ZR_env_norm.data = (R_env_norm.data * Z_env_norm.data) 
	ZR_env_norm.normalize()

	char_func = cross_corr_trace.copy()
	char_func.data = (abs(cross_corr_trace.data) * ZR_env_norm.data) 
	char_func_original = cross_corr_trace.copy()
	char_func_original.data = (cross_corr_trace.data * ZR_env_norm.data)

	################ defining time window ######################################
	t_start = 0.0
	t_end = len(char_func.data)
	i_max = np.argmax(char_func.data)
	char_func_max = char_func_original.data[i_max]
	if max(char_func.data) >= threshold:
		found = 'yes'
		for i in range(i_max,len(char_func.data)):
			if char_func.data[i] <= max(char_func.data) * 0.75:
				i_cc_end = i
				break
		
		t_end = i_cc_end * delta
		for j in reversed(range(0, i_max)):
			if char_func.data[j] <= max(char_func.data) * 0.75:
				i_cc_start = j
				break

		t_start = i_cc_start * delta
			
	else:
		found = 'no'
	

	############## Compute mean ellipticity inside the time window #######################
	Z_R = Z_env.copy()	
	Z_R.data = R_env.data / Z_env.data	
	if found == 'yes':		 
		rayleigh_start = t_start
		rayleigh_start_index = int(rayleigh_start / delta)
		rayleigh_stop = t_end
		rayleigh_stop_index = int(rayleigh_stop / delta)

		Z_R_cut = []
		
		E = 0.0
		sd = 0.0
					
		for i in range(rayleigh_start_index, rayleigh_stop_index):
			Z_R_cut.append(Z_R.data[i])
			
		Z_R_cut_mean = np.mean(Z_R_cut)
		Z_R_cut_sd = np.std(Z_R_cut)

		E = Z_R_cut_mean
		sd = Z_R_cut_sd

		sd_perc = 100 * sd / E
		
					
									
	else:
		rayleigh_start = 0.0
		rayleigh_stop = 0.0
		rayleigh_start_index = 0.0
		rayleigh_stop_index = 0.0
		E = 'NOT FOUND'
		sd = 'NOT FOUND'
	
################################  Plotting  ##############################
	if plot == 'yes':
			if found == 'yes':
				t1 = t_start/60.0 - 50*T0/60.0
				t2 = t_end/60.0 + 50*T0/60.0
			else:
				t1 = 0.0
				t2 = 300.0
			plt.figure(2, figsize=(8.27,11.69))
			plt.rcParams['xtick.labelsize'] = 20
			plt.rcParams['ytick.labelsize'] = 18
			t = 1000
			plt.subplot(511)
			plt.locator_params(axis="y", nbins=5)
			az, dist = take_az_dist(event, station,  cmtpath)
			mw = take_magnitude(event,cmtpath)
			plt.subplots_adjust(left = None, bottom = 0.07, right = None, top = 0.87, wspace = None, hspace = 0.5)
			if found == 'yes':				
				plt.suptitle("Station: "+ station + ' ' + 'Period: ' + str(round(T0,1)) + 's ' + 'Ellipticity: ' + str(round(E,3)) + '\n' + \
						'Event: '+ event.split("_")[1] + "  Mw: " + str(mw) + " Distance: " + str(round(dist,1)) + "deg" , y=0.96, fontsize=20)

			else:
				plt.suptitle(station + ' '  + channel + ' ' + 'period: ' + str(int(Tmin)) + 's - ' + str(int(Tmax)) + 's' +'\n' \
									+ 'Ellipticity: ' + E, y=0.96)
			plt.ticklabel_format(style = 'sci', axis='y', scilimits=(0,0))
			plt.title('Filtered H and V components', fontsize=20 )
			plt.ylabel("Displacement (m)", fontsize=15)
			plt.plot( R_tr.times()/60.0, R_tr.data, linewidth = 2,  color = "black", label = 'H filtered')
			plt.plot( Z_tr.times()/60.0, Z_tr.data, linewidth = 2, color = "red", label = 'V filtered')
			plt.xlim(t1,t2)
			plt.legend(fontsize = 10)

			plt.subplot(512)
			plt.locator_params(axis="y", nbins=5)
			plt.ylabel("Displacement (m)", fontsize=15)
			plt.ticklabel_format(style = 'sci', axis='y', scilimits=(0,0))
			plt.title('Shifted V component', fontsize=20)
			plt.plot( R_tr.times()/60.0, R_tr.data, linewidth = 2,label = 'H filtered',  color = "black")
			plt.plot( Z_dephase.times()/60.0, Z_dephase.data, linewidth = 2, color = "red", label = 'V filtered and shifted')
			plt.xlim(t1,t2)
			plt.legend(fontsize = 10)

			plt.subplot(513)
			plt.locator_params(axis="y", nbins=5)
			plt.title('Cross correlation and normalized HV envelope', fontsize=20)
			plt.plot( cross_corr_trace.times()/60.0, cross_corr_trace.data, linewidth = 2.0, color = "black", label = 'RZ Cross correlation')
			plt.plot( ZR_env_norm.times()/60.0, ZR_env_norm.data, linewidth = 2.0, label = 'Normalized H*V envelope',  color = "black", linestyle = "--")
			plt.xlim(t1,t2)
			plt.legend(fontsize = 10, loc=4)

			plt.subplot(514)	
			plt.locator_params(axis="y", nbins=5)
			plt.title('Characteristic function', fontsize=20)
			#plt.plot(cross_corr_trace.times()/60.0, cross_corr_trace.data )
			plt.plot(char_func.times()/60.0, char_func.data, color = "black", linewidth = 2.0)
			plt.axvline(t_start/60.0, color = "black")
			plt.axvline(t_end/60.0, color = "black")
			plt.xlim(t1,t2)
			plt.ylim(-1,1)
	
			plt.subplot(515)
			plt.title("H/V envelope ratio", fontsize=20)
			plt.plot(Z_R.times()/60.0, Z_R.data, linewidth=2, color = "black")
			plt.ylim(0,5)
			plt.axvline(t_start/60.0, color = "black")
			plt.axvline(t_end/60.0, color = "black")
			plt.xlim(t1,t2)

		#	if found == 'yes':
				#plt.axhline(ell)
				#plt.axhline(Z_R_cut_mean, color = "green")
			#	plt.axhline(Z_R_cut_mean + Z_R_cut_sd, color = 'red')	
			#	plt.axhline(Z_R_cut_mean - Z_R_cut_sd, color = 'red')	
			plt.xlabel('Time (minutes)', fontsize=20)
			plt.locator_params(axis="y", nbins=5)
			plt.savefig(station + '_' + str(start_time) + '_' +  str(T0) + 'k'+str(k)+'_scheme_II.png')
			plt.close()
						
	#except:
	#	print 'Problem! Check: ', event, '\n\n' 
	#	break
	
#	E = np.log10(E)
	
	return T0, E, sd, snr, char_func_max













