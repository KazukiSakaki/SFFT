import numpy as np
import matplotlib.pyplot as plt
import os


def SFFT(T,target,name):
	#experiment parameters
	combustion_time = np.where((T>=0.0)&(T<=3.0))[0]
	steady_time = 1.0

	#FFT parameters
	fs = 10000.0 #sampling rate [Hz]
	dt = 1/fs #time interval [s]
	N = 1024 #number of samples
	df = fs/N #resolution of frequency

	start = min(combustion_time)
	end = max(combustion_time)
	FFT_startTime = 0.0
	FFT_start = np.where((T>FFT_startTime-1/fs)&(T<FFT_startTime+1/fs))[0][0]
	#print FFT_start
	FFT_end = FFT_start + N
	original_FFT_start = FFT_start

	freq_list = np.fft.fftfreq(N,1.0/fs)
	freq_list = freq_list[freq_list>=0]
	#window function (Hamming)
	hammin_window = np.hamming(N)

	#FFT calculation for each time
	target_PSD_result = np.array([])
	steady_PSD_list = np.array([])
	FFT_time = np.array([FFT_startTime])
	time = FFT_startTime
	counter = 0
	steadyflag = 0
	steadyStartCounter = 0
	while  time < 3.0:
		if(steadyflag == 0 and time >3.0):
			steadyflag = 1
			steadyStartCounter = counter
		T_FFT = T[FFT_start:FFT_end]
		target_FFT = target[FFT_start:FFT_end]
		target_FFT = target_FFT - np.average(target_FFT) #subtract average value (to remove DC components)
		#complex result
		target_FFT_result = np.fft.fft(target_FFT*hammin_window)
		#Power Spectral Density of Pc
		target_PSD = np.power(np.absolute(target_FFT_result),2)/df
		if time == FFT_startTime:
			frequencyMap = target_PSD[0:N/2]
		if time > FFT_startTime:
			steady_PSD_list = np.append(steady_PSD_list,target_PSD)

		frequencyMap = np.vstack((frequencyMap,target_PSD[0:N/2]))

		#time advance
		FFT_start = FFT_start + N
		FFT_end = FFT_start + N
		time = (FFT_start- original_FFT_start)*dt + FFT_startTime
		FFT_time = np.append(FFT_time,time)

		counter = counter + 1

	#Graph drawing
	frequencyMap = frequencyMap.T
	PSDmax = max(steady_PSD_list[N/2*steadyStartCounter:])
	PSDmin = min(steady_PSD_list)
	interval = np.linspace(0,PSDmax,1000)
	plt.figure()
	CF = plt.contourf(FFT_time,freq_list,frequencyMap,interval,cmap=plt.cm.jet)
	CB = plt.colorbar(CF)
	CB.set_label("Power Spectrum Density[MPa$^2$/Hz]")
	#plt.title('Spectrogram of '+testname)
	plt.xlabel('Time [s]')
	plt.ylabel('Frequency [Hz]')
	plt.savefig(testname + "_" +name+ "_FFT"+".png")
	#plt.show()
	plt.close()

	#logarithmic version
	frequencyMap_dB = 10*np.log10(frequencyMap)
	PSDmax_dB = 10*np.log10(PSDmax)
	PSDmin_dB = 10*np.log10(PSDmin)
	interval = np.linspace(PSDmin_dB,PSDmax_dB,1000)
	plt.figure()
	CF = plt.contourf(FFT_time,freq_list,frequencyMap_dB,interval,cmap=plt.cm.jet)
	CB = plt.colorbar(CF)
	CB.set_label("Power Spectrum Density [dB]")
	#plt.title('Spectrogram of '+testname)
	plt.xlabel('Time [s]')
	plt.ylabel('Frequency [Hz]')
	plt.savefig(testname +"_" +name+ "_FFT_log"+".png")
	#plt.show()
	plt.close()

runnumber = 38
stop = 45
while runnumber < stop:
	testname = "RawData/AX-OC-RUN" + str(runnumber)
	extention = ".CSV"

	if os.path.isfile(testname+extention):
		print testname +" is under calculation"
		data = np.genfromtxt(testname+extention,delimiter=",",skip_header=14)

		T=data[:,0]-1.0
		F = data[:,1]
		PTP = data[:,2]
		PFT = data[:,3]
		POT = data[:,4]
		PFI = data[:,5]
		POI = data[:,6]
		PFL = data[:,7]
		POL = data[:,8]
		PFFL = data[:,9]
		POFL = data[:,10]
		PcS = data[:,11]
		#PcD = Pc
		PcD = data[:,12]
		TRG = data[:,13]
		TIG = data[:,14]
		TFI = data[:,15]
		TOI = data[:,16]

		SFFT(T,PcD,"PcD")
		SFFT(T,PFI,"PFI")
		SFFT(T,POI,"POI")
		runnumber = runnumber + 1
	else:
		print testname +" does not exist"
		runnumber = runnumber + 1


def newFunction():
	return
	
