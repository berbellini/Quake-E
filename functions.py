#############################################
#         				    				#
#   Andrea Berbellini - PhD student         #
#       INGV-Bologna   September 2019 	    #
#	andrea.berbellini@ingv.it           	#
#					    					#
#############################################

import numpy as np
from obspy.geodetics.base import gps2dist_azimuth
import scipy
import math



def sig_noise(tr, s):
	tr_c = tr.copy()
	delta = tr_c.stats.delta
	tr_c.data = abs(tr_c.data)
	mean_noise = np.mean(tr_c.data[0: int(s / delta)])
	max_signal = max(tr_c.data)
	sig_noise = max_signal / mean_noise

	return sig_noise

#-------------------------------------------------------------------

def phase_shift(trace, shift):
	trace_shift = trace.copy()
	trace_long = trace.copy()
	shift_points = int(shift / trace_shift.stats.delta)
	if shift >= 0.0:
		for i in range(0, len(trace.data) - shift_points ):
			trace_shift.data[i] = trace.data[i + shift_points]
	else:
		
		trace_shift = add_zero(trace,abs(shift) )
		trace_shift.trim(trace.stats.starttime, trace.stats.endtime )
	
	return trace_shift

#--------------------------------------------------------------------

def cross_correlation(trace1, trace2, T_window):
	T_window = T_window / trace1.stats.delta
	s1 = trace1.copy()
	s1_quad = s1.copy()
	s1_quad.data = np.square(s1.data) 
	s1_quad_integral = s1_quad.copy()
	s1_quad_integral = scipy.integrate.cumtrapz(s1_quad_integral)

	s2 = trace2.copy()
	s2_quad = s2.copy()
	s2_quad.data = np.square(s2.data) 
	s2_quad_integral = s2_quad.copy()
	s2_quad_integral = scipy.integrate.cumtrapz(s2_quad_integral)

	minsamples = min(len(s1), len(s2))
	a1a2 = s1.copy()
	for i in range(0,min(len(s1),len(s2))):
		a1a2.data[i] = s1.data[i] * s2.data[i]

	

	A = np.arange(0,len(trace1)-T_window, T_window)
	a_med = []	
	cross_corr = []

	for a in A:
		t1 = a 
		t2 = a + T_window
		a1a2_integral = integral(a1a2, t1, t2)
		s1_integral = integral(s1_quad, t1, t2)
		s2_integral = integral(s2_quad, t1, t2)
		time = (a + (T_window/2.0) ) * trace1.stats.delta
		a_med.append(time)
		corr =  a1a2_integral/( math.sqrt(s1_integral) * (math.sqrt(s2_integral)) )
		cross_corr.append(corr)
	
	cross_corr_trace = trace1.copy()

	cross_corr_trace.data = np.interp( cross_corr_trace.times(), a_med, cross_corr )	

	
	return cross_corr_trace

#---------------------------------------------------------------------------------------
def integral(trace, t1, t2):
	a = trace.copy()
	interval = np.arange(t1,t2, a.stats.delta)
	a_int = 0.0
	for i in interval:
		a_int = a_int + (a.data[int(i)] * a.stats.delta)
	return a_int


#---------------------------------------------------------------------------------------


def take_station_coordinate(station):
	#station_file = "/data/berbellini/Copy/Ellipticity_UCL/lib/stations.txt"
	station_file = "./stations.txt"
	#print station_file, station
	lines = open(station_file, "r").readlines()
	for i in range(2,len(lines)):
		if lines[i].split()[0] == station:
			lat = float(lines[i].split()[2])
			lon = float(lines[i].split()[3])
	return lat, lon



def take_az_dist(event, station,  cmtpath):
	R = 6371.0
	cmtfile_lines = open(cmtpath + "/" + event, "r").readlines()
	ev_lat = float(cmtfile_lines[4].split()[1])
	ev_lon = float(cmtfile_lines[5].split()[1])
	stat_lat, stat_lon = take_station_coordinate(station)
	dist = gps2dist_azimuth(stat_lat, stat_lon, ev_lat, ev_lon)[0] / 1000.0
	baz = gps2dist_azimuth(stat_lat, stat_lon, ev_lat, ev_lon)[1]
	az = gps2dist_azimuth(stat_lat, stat_lon, ev_lat, ev_lon)[2]
	dist_grad = (dist / R) * 360 / (2*np.pi) 
	return baz, dist_grad


def take_magnitude(event,cmtpath):
	filename = cmtpath + '/' + event
	fl = open(filename, "r")
	lines = fl.readlines()
	if lines[0].split()[0] == "PDE":
		mw = float(lines[0].split()[11])
	elif lines[0].split()[0][0:4] == "PDEW":
		mw = float(lines[0].split()[10])
	return mw











