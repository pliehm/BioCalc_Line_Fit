#cython: boundscheck=False
import numpy as np
cimport numpy as np
cimport cython


DTYPE = np.uint16
ctypedef np.uint16_t DTYPE_t



# fit simulated and experimental minima to calculate thickness
# problem: speed!! --> one main problem is the access of the large list sim_waves --> can this be written as an array?


@cython.boundscheck(False)
cdef Fit(list sim_waves, list exp_waves, int tolerance):

# cython type definition of the used objects
 
	cdef list sim_min_waves = [[1000],[1000]] # dummy values to have two lists
	cdef unsigned short i, k
	cdef unsigned short L_exp_waves = len(exp_waves) # too not use the len() function too often
	cdef  float summe=0
	cdef list results=[]

	if L_exp_waves > 2: # first test if there are more than one minima
		for i in xrange(len(sim_waves)): # do the following calculations for every simulated thickness
			if sim_waves[i][1] == L_exp_waves: # case for equal minimas (exp, sim)
				summe=0
# perform something like least-square with every exp-wavelength 
				for k in xrange(len(exp_waves)): 
					summe+= abs(sim_waves[i][2][k]-exp_waves[k])
# append the thickness and error to a list
				sim_min_waves[0].append(sim_waves[i][0])
				sim_min_waves[1].append(summe/float(len(exp_waves)))
# do the same if number of exp and sim is not equal
			if sim_waves[i][1] == (L_exp_waves + 1):
				summe=0
# check if the first elements (exp and sim) or the last tow are not matching
				if abs(sim_waves[i][2][0] - exp_waves[0]) > abs(sim_waves[i][2][-1]-exp_waves[-1]):
					for k in xrange(len(exp_waves)):
						summe+= abs(sim_waves[i][2][k+1]-exp_waves[k])
					sim_min_waves[0].append(sim_waves[i][0])
					sim_min_waves[1].append(summe/float(len(exp_waves)))
				else:
					for k in xrange(len(exp_waves)):
						summe+= abs(sim_waves[i][2][k]-exp_waves[k])
					sim_min_waves[0].append(sim_waves[i][0])
					sim_min_waves[1].append(summe/float(len(exp_waves)))
			if sim_waves[i][1] == (L_exp_waves - 1):

				summe=0
				if abs(sim_waves[i][2][0] - exp_waves[0]) > abs(sim_waves[i][2][-1]-exp_waves[-1]):
					for k in xrange(sim_waves[i][1]):
						summe+= abs(sim_waves[i][2][k]-exp_waves[k+1])
					sim_min_waves[0].append(sim_waves[i][0])
					sim_min_waves[1].append(summe/float(len(exp_waves)))
				else:
					for k in xrange(sim_waves[i][1]):
						summe+= abs(sim_waves[i][2][k]-exp_waves[k])
					sim_min_waves[0].append(sim_waves[i][0])
					sim_min_waves[1].append(summe/float(len(exp_waves)))
				
# append results to a list for further processing
		#results=[]
		#results=[sim_min_waves[i][1] for  i in xrange(len(sim_min_waves))]

# return the thickness with minimum value
#		if  len(sim_min_waves[1])>1 and (min(results) < tolerance):
#			return sim_min_waves[ results.index(min(results))][0]
		if  len(sim_min_waves[0])>1 and (min(sim_min_waves[1]) < tolerance):
			return sim_min_waves[0][sim_min_waves[1].index(min(sim_min_waves[1]))]
		else: 
			return 0

	else:
		return 0

# function to get minima of one array


cdef peakdetect(y_axis, x_axis = None, unsigned short lookahead_min=5, unsigned short lookahead_max=3, unsigned short delta = 0):
	
	# define output container
	#cdef list max_peaks=[]
	cdef list min_peaks = []
	cdef list dump = [] # used to pop the first hit which almost always is false

	# check input data --> this makes the algorithm 5 times slower
	#x_axis, y_axis = _datacheck_peakdetect(x_axis, y_axis) 
	
	# store data length for later use
	cdef int length = len(y_axis)

	#perform some checks
	#if lookahead < 1:
	#	raise ValueError, "Lookahead must be '1' or above in value"
	#if not (np.isscalar(delta) and delta >= 0):
	#	raise ValueError, "delta must be a positive number"

	#maxima and minima candidates are temporarily stored in 
	#mx and mn respectively
	cdef int mn = 1000
	cdef int mx = -1000
	cdef unsigned short x,y, index

	for index, (x,y) in enumerate(zip(x_axis[:-lookahead_min], y_axis[:-lookahead_min])):
		
		if y > mx:
			mx = y
			mxpos = x

		if y < mn:
			mn = y
			mnpos = x

		#### look for max ####
		
		if y < mx-delta and mx != 1000:
			#Maxima peak candidate found
			# lool ahead in signal to ensure that this is a peak and not jitter
			if y_axis[index:index+lookahead_max].max() < mx:
				#max_peaks.append([mxpos, mx])
				dump.append(True)
				#set algorithm to only find minima now
				
				mx = 1000
				mn = 1000
				if index+lookahead_min >= length:
					#end is within lookahead no more peaks can be found
					break
				continue

		#### look for min ####	
		
		if y > mn+delta and mn != -1000:
			#Minima peak candidate found
			# look ahead in signal to ensure that this is a peak and not jitter
			if y_axis[index:index+lookahead_min].min() > mn:
				min_peaks.append(mnpos)
				dump.append(False)
				#set algorithm to only find maximum now
				mn = -1000
				mx = -1000
				if index+lookahead_min >= length:
					#end is within lookahead no more peaks can be found
					break


	#Remove the false hit on the first value of the y_axis
	if len(dump)>0:
		if not dump[0]:
			min_peaks.pop(0)
		else:
			pass	
	
		#no peaks were found, should the function return empty lists?
	
	return min_peaks

	
# find reflection minima for every pixel

def c_Fit_Pixel(np.ndarray[DTYPE_t, ndim=3] data, int zeile, list sim_waves, list waves, int tolerance, unsigned short lookahead_min,unsigned short lookahead_max, unsigned short delta):
	cdef unsigned short spalte
	cdef np.ndarray[DTYPE_t, ndim=1] intensity
	cdef list minima_exp=[]
	cdef list thickness_ready=[]
	cdef int counter=0
	for spalte in xrange(1280):
		intensity = data[:,zeile, spalte]
		minima_exp = peakdetect(intensity, waves, lookahead_min,lookahead_max, delta)
		#print minima_exp
		thickness_ready.append(Fit(sim_waves, minima_exp,tolerance))
		#print counter		
		counter+=1
		
	return thickness_ready
	

