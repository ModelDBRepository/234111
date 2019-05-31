#!/opt/local/bin/python2.7
# -*- coding: utf-8 -*-

import sys, getopt, os
import code  #<---cause script to go interactive at any line by inserting code.interact(local=locals())
import time

import neuron
import nrn

from neuron import *
import matplotlib
matplotlib.rcParams['agg.path.chunksize'] = 10000

import matplotlib.pyplot as pyplot
import pickle
import numpy as np

from neuron import h
h.load_file('stdrun.hoc')    # defines h.cvode
home = os.getenv("HOME")
current_dir = os.getcwd()+'/'
folder = current_dir



##############################################################################
################# Initialize the Model #######################################

###Specify directory to load user-defined mechanisms (e.g. CLS.mod) 
# load_mechanisms(os.getenv("HOME")+'')


### Create a section (node of ranvier):
ranvier = h.Section(name='ranvier')
h.psection(sec=ranvier)
# print dir(ranvier)

###All in dimensions of micrometre^n
VolumeOut = 3.0

Area = 6.0
VolumeIn = 3.0
radius = 2.0*VolumeIn/Area
length = (Area**2.0)/(4.0*np.pi*VolumeIn)
diameter = 2.0*radius

ranvier.L = length
ranvier.diam = diameter

print "diameter of section = ", ranvier.diam
print "Length of section = ", ranvier.L



shape_window = h.PlotShape()
shape_window.exec_menu('Show Diam')

### Biophysics
for sec in h.allsec():
    sec.Ra = 100    # Axial resistance in Ohm * cm
    sec.cm = 1      # Membrane capacitance in micro Farads / cm^2

### Insert active CLS Hodgkin-Huxley current with pumps in the ranvier
ranvier.insert('CLS')
cvode = h.CVode()
cvode.active(1)
cvode.use_long_double(1)
cvode.atol(1.0e-6)
# h.finitialize(v0)



ranvier.AC_CLS = 0.0 ###Proportion of affected (left-shifted) SODIUM channels on node
ranvier.ACpotassium_CLS = 0.0 ###Proportion of affected (left-shifted) POTASSIUM channels on node
ranvier.vLeftShift_CLS = 0.0			#Coupled Left-Shift voltage (mV)
ranvier.timeLS_CLS = 220.0			#Delay (ms) prior to activation of coupled Left-Shift 
print  "Voltage Left-Shift = ", ranvier.vLeftShift_CLS


# ranvier.gnabar_CLS = 0.120			# Maximum specific sodium channel conductance in S/cm2
# ranvier.gkbar_CLS = 0.036			# Maximum potassium channel conductance in S/cm2
ranvier.ena = 51.5				# The reversal potential for the sodium channel [Default value = 50 mV]
ranvier.ek = -81.3					# The reversal potential for the potassium channel [Default value = -77 mV]
ranvier.gl_CLS = 0.0005 			# Leakage conductance in S/cm2
ranvier.el_CLS = -59.9 #-51.5		# Leakage reversal potential in mV
# ranvier.m_CLS                 	# sodium activation state variable
# ranvier.h_CLS                		# sodium inactivation state variable
# ranvier.n_CLS                 	# potassium activation state variable
# ranvier.ina_CLS         			# sodium current through the hh channels in mA/cm2
# ranvier.ik_CLS       	  			# potassium current through the hh channels mA/cm2


h.celsius = 20.0
ranvier.nai = 20.0
ranvier.nao = 154.0
ranvier.ki = 150.0
ranvier.ko = 6.0
ranvier.v = -6.00724451e+01  #-59.8137886



h.psection(sec=ranvier)
stim = h.IClamp(ranvier(0.5))
stim.delay = 200.0*1e03
stim.dur = 20.0*1e03


##############################################################################
################# Run the simulation #########################################

def runAC_vLS(AC=1.0, vLS=2.0, Istim=12.0, duration=300.0, OpenTheFigures=False, LoadFromPrevious=False, SaveData=False, ClipLeft=3):



	#Istim = 12.0 (microamps per cm^2)
	stim.amp = (Istim*1.0E-05)*Area

	ranvier.AC_CLS = AC ###Proportion of affected (left-shifted) SODIUM channels on node
	ranvier.ACpotassium_CLS = 0.0 ###Proportion of affected (left-shifted) POTASSIUM channels on node
	ranvier.vLeftShift_CLS = vLS			#Coupled Left-Shift voltage (mV)


	filename = 'AC='+str(ranvier.AC_CLS)+'vLS='+str(ranvier.vLeftShift_CLS)+'Istim='+str(Istim)

	
	if LoadFromPrevious == True:

		try:
			# data = np.load(folder+filename+'.npy')
			# data = data.item() #<---Recover dictionary

			data = np.load(folder+filename+'.npz')
			data = data.items()[0][1].item() #<---Recover dictionary

			t_vec = data['t_vec']
			v_vec = data['v_vec']
			Istim = data['Istim']
			StimDelay = data['StimDelay']
			StimDuration = data['StimDuration']


			fi_data, spikes, freqVsTime = spike_detection(v_vec, t_vec, Istim, StimDelay=StimDelay, StimDuration=StimDuration, ClipLeft=ClipLeft) #<---get frequency and spike data versus current
			data['fi_data'] = fi_data
			data['spikes'] = spikes
			data['freqVsTime'] = freqVsTime
		except:
			print 'Failed either to locate or load data: running the simulation.'
			LoadFromPrevious==False


	if LoadFromPrevious==False:

		simdur = duration*1e03
		h.tstop = simdur

		### Reset up the recording vectors, run, and plot.
		t_vec = h.Vector()        # Time stamp vector
		i_vec = h.Vector()
		vLS_vec = h.Vector()        # Membrane potential vector
		v_vec = h.Vector()  
		ena_vec = h.Vector() 
		ek_vec = h.Vector()
		nao_vec = h.Vector() 
		nai_vec = h.Vector()
		ko_vec = h.Vector() 
		ki_vec = h.Vector()
		ina_vec = h.Vector() 
		ik_vec = h.Vector()
		ink_vec = h.Vector()  #<---pump
		m_vec = h.Vector() 
		h_vec = h.Vector() 
		n_vec = h.Vector() 
		mLS_vec = h.Vector() 
		hLS_vec = h.Vector()  
		nLS_vec = h.Vector() 
		 


		i_vec.record(stim._ref_i)
		# t_vec.record(h._ref_t)
		# v_vec.record(ranvier(0.5)._ref_v)
		cvode.record(ranvier(0.5)._ref_v, v_vec, t_vec)
		vLS_vec.record(ranvier(0.5)._ref_vLS_CLS)
		ena_vec.record(ranvier(0.5)._ref_ena)
		ek_vec.record(ranvier(0.5)._ref_ek)
		nao_vec.record(ranvier(0.5)._ref_nao)
		nai_vec.record(ranvier(0.5)._ref_nai)
		ko_vec.record(ranvier(0.5)._ref_ko)
		ki_vec.record(ranvier(0.5)._ref_ki)
		ina_vec.record(ranvier(0.5)._ref_ina)
		ik_vec.record(ranvier(0.5)._ref_ik)
		ink_vec.record(ranvier(0.5)._ref_ink_CLS)
		m_vec.record(ranvier(0.5)._ref_m_CLS)
		h_vec.record(ranvier(0.5)._ref_h_CLS)
		n_vec.record(ranvier(0.5)._ref_n_CLS)
		mLS_vec.record(ranvier(0.5)._ref_mLS_CLS)
		hLS_vec.record(ranvier(0.5)._ref_hLS_CLS)
		nLS_vec.record(ranvier(0.5)._ref_nLS_CLS)
		# m_vec, h_vec, n_vec	

		# h.steps_per_ms = 1000
		# h.dt = 1.0/h.steps_per_ms

		#**************************************
		start = time.time()#*******************

		h.run()
		# stim.delay = 1050.0*1e03
		# h.continuerun(1200*1e3)

		end = time.time()#*********************
		#**************************************



		RunTime = end - start
		print '\n\n\n' + 'RunTime = ', RunTime, ' seconds \n\n\n'

		##############################################################################
		################# repackage the data: convert HOC vectors to numpy arrays ####

		t_vec = np.array(t_vec)
		i_vec = np.array(i_vec)
		v_vec = np.array(v_vec)
		vLS_vec = np.array(vLS_vec)
		ena_vec = np.array(ena_vec)
		ek_vec = np.array(ek_vec)
		nao_vec = np.array(nao_vec)
		nai_vec = np.array(nai_vec)
		ko_vec = np.array(ko_vec)
		ki_vec = np.array(ki_vec)
		ina_vec = np.array(ina_vec)
		ik_vec = np.array(ik_vec)
		ink_vec = np.array(ink_vec)
		m_vec = np.array(m_vec)
		h_vec = np.array(h_vec)
		n_vec = np.array(n_vec)
		mLS_vec = np.array(mLS_vec)
		hLS_vec = np.array(hLS_vec)
		nLS_vec = np.array(nLS_vec)
		
		timeSeries = np.array([t_vec, i_vec, v_vec, vLS_vec, ena_vec, ek_vec, nao_vec, nai_vec, ko_vec, ki_vec, ina_vec, ik_vec, ink_vec, m_vec, h_vec, n_vec, mLS_vec, hLS_vec, nLS_vec])
		names = np.array(['t_vec', 'i_vec', 'v_vec', 'vLS_vec', 'ena_vec', 'ek_vec', 'nao_vec', 'nai_vec', 'ko_vec', 'ki_vec', 'ina_vec', 'ik_vec', 'ink_vec', 'm_vec', 'h_vec', 'n_vec', 'mLS_vec', 'hLS_vec', 'nLS_vec'])

		fi_data, spikes, freqVsTime = spike_detection(v_vec, t_vec, Istim, StimDelay=stim.delay, StimDuration=stim.dur, ClipLeft=ClipLeft) #<---get frequency and spike data versus current


		data = dict(Istim = Istim, StimDelay=stim.delay, StimDuration=stim.dur, AC = ranvier.AC_CLS, vLS = ranvier.vLeftShift_CLS, t_vec = t_vec, i_vec=i_vec, v_vec = v_vec, vLS_vec = vLS_vec, ena_vec = ena_vec, ek_vec = ek_vec, nao_vec = nao_vec, nai_vec = nai_vec, ko_vec = ko_vec, ki_vec = ki_vec, ina_vec = ina_vec, ik_vec = ik_vec, ink_vec = ink_vec, m_vec = m_vec, h_vec = h_vec, n_vec = n_vec, mLS_vec = mLS_vec, hLS_vec = hLS_vec, nLS_vec = nLS_vec, fi_data = fi_data, spikes = spikes, freqVsTime = freqVsTime)

		###Save the data
		if SaveData == True:
			# np.save(folder+filename, data)
			np.savez_compressed(folder+filename, data)

		if OpenTheFigures == True:
			plot(data)



		### Print initial and final state of system
		for i in np.arange(len(timeSeries)):
			print names[i], ' = ' , timeSeries[i,1]
			# print names[i], (timeSeries[i])[-1]
		print '\n\n'

		for i in np.arange(len(timeSeries)):
			print names[i], ' = ' , timeSeries[i,-1]



		print 'stim.delay = ', stim.delay
		print 'stim.dur = ', stim.dur


	return data




##############################################################################
################# Frequency Calculation ######################################

from scipy.signal import argrelextrema


def spike_detection(Voltage, Time, Injected_Current, StimDelay, StimDuration, ClipLeft=1, MinAmplitude=5.0):


	# for local maxima
	locMax = argrelextrema(Voltage, np.greater)[0]

	# for local minima
	locMin = argrelextrema(Voltage, np.less)[0]

	print 'len(locMax), len(locMin) = ', len(locMax), len(locMin)
	minlen = min(len(locMax), len(locMin))
	locMax, locMin = locMax[:minlen], locMin[:minlen] 


	
	AmplitudeCutoff = np.where((Voltage[locMax] - Voltage[locMin]) > MinAmplitude)
	locMax = locMax[AmplitudeCutoff]
	locMin = locMin[AmplitudeCutoff]


	Vspikes = np.array(Voltage[locMax])
	Tspikes = np.array(Time[locMax])

	TimeCutoff = np.where((Tspikes > StimDelay) & (Tspikes < StimDelay+20.0*1e03))
	Vspikes = Vspikes[TimeCutoff][ClipLeft:]
	Tspikes = Tspikes[TimeCutoff][ClipLeft:]

	deltaTspikes = (Tspikes[1:] - Tspikes[:-1])*1e-3 #<---Exponent converts calcultated frequency to Hz

	# select = np.where((Tspikes[1:] - Tspikes[:-1]) < 100.0)
	# Vspikes = Vspikes[select]
	# Tspikes = Tspikes[select]


	runnningFreq = deltaTspikes**(-1.0)
	FreqCutoff = np.where(runnningFreq > 20.0)
	freqVsTime = np.array((runnningFreq[FreqCutoff], Tspikes[1:][FreqCutoff]))


	spikes = np.array((Vspikes, Tspikes))


	N = np.array(len(spikes[0])).astype('float64') #<--- Number of spikes counted
	print 'N = ', N
	if N > 1:
		Spiking_Detected = True
		t1 = spikes[1,0]
		t2 = spikes[1,-1]
		# t1 = Time[spikes[1,0].astype(int)] #<--- Time of first (counted) spike, in milliseconds
		# t2 = Time[spikes[1,-1].astype(int)] #<--- Time of last (counted) spike, in milliseconds

		frequency = N/(t2-t1)
		frequency *= 1.0e3 #<---Convert from millihertz to hertz
	else:
		Spiking_Detected = False
		print 'WARNING: NO SPIKING DETECTED'
		frequency = 0.0
		
	print 'frequency at ', Injected_Current, ' microamps per cm2 ' , ' = ', frequency, 'Hz'
	return np.array([Injected_Current, frequency]), spikes, freqVsTime





#### Program can be used as either a SCRIPT OR MODULE, since code below will activate only when program is run from the terminal
#### and this code does nothing when program is imported as a module:
if __name__ == "__main__":
    import sys
    import numpy as np
    # print 'sys.argv[0] = ', sys.argv[0]
	# print 'sys.argv[1] = ', sys.argv[1]
	# print 'sys.argv[2] = ', sys.argv[2]

    data_dict = bbm.runAC_vLS(AC=sys.argv[1], vLS=sys.argv[2], duration=sys.argv[3])






intit = runAC_vLS(duration=1.0, SaveData=False)
print '\n\n\n' + '******** INITIALIZATION COMPLETE ********' + '\n\n\n'


##############################################################################
################# Plotting ###################################################
def plot(data, xmin=-10.0, xmax=None, skip=1, OpenTheFigures=True, InteractiveFigure=False):

	fig, (ax1, ax2) = pyplot.subplots(nrows=2, ncols=1, figsize=(8,8), sharex=True)
	figurename = 'AC='+str(np.round(data['AC'], 3))+'vLS='+str(np.round(data['vLS'], 3))+'Istim='+str(np.round(data['Istim'], 3))
	ax1.set_title(figurename)

	
	colour = 'magenta' 


	mLS_vec  = data['mLS_vec']
	StimDuration  = data['StimDuration']
	freqVsTime  = data['freqVsTime']
	ink_vec  = data['ink_vec']
	AC  = data['AC']
	ina_vec  = data['ina_vec']
	Istim  = data['Istim']
	v_vec  = data['v_vec']
	vLS_vec  = data['vLS_vec']
	ik_vec  = data['ik_vec']
	fi_data  = data['fi_data']
	vLS  = data['vLS']
	n_vec  = data['n_vec']
	nao_vec  = data['nao_vec']
	ko_vec  = data['ko_vec']
	nai_vec  = data['nai_vec']
	ki_vec  = data['ki_vec']
	h_vec  = data['h_vec']
	spikes  = data['spikes']
	i_vec  = data['i_vec']
	hLS_vec  = data['hLS_vec']
	nLS_vec  = data['nLS_vec']
	t_vec  = data['t_vec']
	StimDelay  = data['StimDelay']
	ena_vec  = data['ena_vec']
	m_vec  = data['m_vec']
	ek_vec  = data['ek_vec']

	

	ax1.plot(t_vec[::skip]*1e-3, ena_vec[::skip], label=r'$E_{Na}$', linewidth=1.0, marker='None')
	ax1.plot(spikes[1]*1e-3, spikes[0], label=r'peaks', linewidth=0, ms=2.0, marker='o', color=colour, mec=colour, mfc=colour, alpha=1.0, zorder=4)
	ax1.plot(t_vec[::skip]*1e-3, v_vec[::skip], label=r'$V_{m}$', linewidth=0.2, marker='None', alpha=0.5)
	# # ax1.plot(t_vec[::skip]*1e-3, vLS_vec[::skip], label='vLS', linewidth=1.0, marker='None')
	ax1.plot(t_vec[::skip]*1e-3, ek_vec[::skip], label=r'$E_{K}$', linewidth=1.0, marker='None')

	# ax1.plot(t_vec[::skip]*1e-3, nao_vec[::skip], label=r'$[Na]_{o}$')
	# ax1.plot(t_vec[::skip]*1e-3, nai_vec[::skip], label=r'$[Na]_{i}$')
	# ax1.plot(t_vec[::skip]*1e-3, ko_vec[::skip], label=r'$[K]_{o}$')
	# ax1.plot(t_vec[::skip]*1e-3, ki_vec[::skip], label=r'$[K]_{i}$')
	# ax1.plot(t_vec[::skip]*1e-3, ina_vec[::skip], label=r'$I_{Na}$', linewidth=0.5, color='orange', alpha=1.0)
	# ax1.plot(t_vec[::skip]*1e-3, ik_vec[::skip], label=r'$I_{K}$', linewidth=0.5, color='blue', alpha=0.6)
	# ax1.plot(t_vec[::skip]*1e-3, ink_vec[::skip], label=r'$I_{pump}$')


	ax2.plot(t_vec*1e-3, i_vec, label=r'$I_{Clamp}$', linewidth=1.0, marker='None', color='black', alpha = 0.6, zorder=4)
	# ax2.plot(freqVsTime[1]*1e-3, freqVsTime[0], label=r'instantaneous frequency', linewidth=0, ms=2.0, marker='o', color=colour, mec=colour, mfc=colour, alpha=1.0, zorder=4)
	


	# ax1.set_xlim([799.9, 810.0])
	ax1.set_xlim([xmin, xmax])
	# ax2.set_ylim(1.1*np.array([min(i_vec), max(i_vec)]))
	# ax2.set_ylim([75.0, 110.0])
	ax1.set_ylabel(r'potential $(mV)$')
	# ax2.set_ylabel(r'frequency $(Hz)$')
	ax2.set_xlabel(r'time $(s)$')
	# ax1.set_ylabel('mV')
	ax1.legend(loc='upper right', fancybox=True, framealpha=0.5)
	ax2.legend(loc='upper right', fancybox=True, framealpha=0.5)


	##### Save and open the figure
	figurePath = folder+figurename
	# pyplot.savefig(figurePath+'.pdf', transparent=True)
	pyplot.savefig(figurePath+'.png', dpi=1000, transparent=False)
	if OpenTheFigures:
		os.system("open " + figurePath+'.png')
		if InteractiveFigure==True:
			pyplot.show()
	pyplot.clf()
	pyplot.close(fig)


	def analyze(data):
		x = 0.0
		return x



