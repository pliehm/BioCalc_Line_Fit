import cython_peakdetect as Fit
import numpy as np
import time
import os 
import Image as im
import matplotlib.pyplot as plt
# enter row you want to readout

zeile = 1022

# enter folder with data

folder = '52_Feature' 

# enter name of simulation_file

sim_file = 'Sim_0.5Cr_20Au_Elastomer_RT601_15Au_500_750nm.txt'

# enter average deviation of experiment to simulation
tolerance=1

# define parameters for peakdetection  
lookahead_min = 3 # something like peak width for thr minima
lookahead_max = 2 # for the maxima --> should not be larger than lookahead_min
delta = 3	# something like peak height

# chose wavelength range and step-width

wave_start = 550
wave_end = 750
wave_step = 1

# chose elastomer thickness range

d_min=9500
d_max=13500

# make wavelength list



waves=[]

waves=[wave_start + i*wave_step for i in xrange((wave_end-wave_start)/wave_step+1)]

## read image data 
dateien=os.listdir(folder)
dateien.sort()

alle=np.zeros((wave_end-wave_start+1,1024,1280),np.uint16)
def image2array(Img):
	newArr= np.fromstring(Img.tostring(),np.uint8)
	newArr= np.reshape(newArr, (1024,1280))
	return newArr

counter=0
for i in xrange(len(dateien)):
	if dateien[i][-5:]=='.tiff':
		if int(dateien[i][:3]) >= wave_start and int(dateien[i][:3]) <= wave_end:
			print dateien[i]
			print counter
			Img=im.open(folder + '/' + dateien[i])
			alle[counter]=image2array(Img)
			counter+= 1


# read simulation file

p= open(sim_file,'r')


string=p.read()

p.close()

positioncounter = 0

sim_waves=[]
wave_block=[]

for thisline in string.split('\n'):
	if ('\t' in thisline) == False and len(thisline) != 0:
		thickness = thisline
	if ('\t' in thisline) == True and int(thickness) >= d_min and int(thickness) <= d_max:
		positioncounter == 0
		for word in thisline.split('\t'):
			if len(word)<6 and float(word) >= wave_start +lookahead_min and float(word)<= wave_end - lookahead_min:
				wave_block.append(float(word))
	if len(thisline) == 0 and int(thickness) >= d_min and int(thickness) <= d_max:
		sim_waves.append([thickness,wave_block])
		wave_block=[]

t1 = time.time()

dicke = Fit.c_Fit_Pixel(alle, zeile, sim_waves, waves, tolerance, lookahead_min, lookahead_max, delta)

t2 = time.time()

print t2-t1, 'Sekunden'

### plot data ###

p = open('boarder_line_fit_row_' + folder +  '_' + str(zeile) + '_tolerance_' + str(tolerance) + '_'+str(wave_start) + '_' + str(wave_end)  + 'nm.txt','w')
for i in range(len(dicke)):
	p.write(str(dicke[i])+'\n')
p.close()	

# convert strings in "dicke" to integers to plot in histogram
# make list without zeros to get real mean value

dicke_i = []

for element in dicke:
	dicke_i.append(int(element))

print 'Anzahl von Dicke = 0', dicke.count(0)


plt.figure(1)
plt.plot(range(len(dicke)),dicke)
plt.axis([0,1280,0,12000])
plt.savefig('both_boarders_' + folder + '_zeile_' +  str(zeile) + '_tolerance_' + str
(tolerance) + '_'+str(wave_start) + '_' + str(wave_end)  + 'nm.png')

plt.figure(2)
plt.hist(dicke_i, bins = 100, color = 'g')
plt.grid()
plt.savefig('hist_' + 'both_boarders_' + folder + '_zeile_' +  str(zeile) + '_tolerance_' + str
(tolerance) + '_'+str(wave_start) + '_' + str(wave_end)  + 'nm.png')

plt.show()



