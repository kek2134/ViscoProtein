import numpy as np
import matplotlib.pyplot as plt
import argparse
import numpy.random as rng
import scipy.interpolate as interp
import numpy.fft as fft




def ComputeCt(x_t):
	dx_t = x_t - np.mean(x_t)
	x_w = fft.fft(dx_t)
	J_w = x_w * np.conjugate(x_w)
	C_t = fft.ifft(J_w)
	return C_t

def ImportXvg(filename, column=[0,1],type=float, max_values=10000):
	dat = []
	frame = 0
	with open(filename) as f:
		for l in f:
			frame += 1
			l_arr = np.array(l.split())
			#print l_arr[column]
			dat.append(np.array([type(val) for val in l_arr[column]]))
			if frame >= max_values:
				break
	return np.array(dat).T

def main():
	data_tr = ImportXvg('./data/dist.xvg', max_values=1E7)
	print data_tr.shape
	F_rt = interp.interp1d(data_tr[0,:], data_tr[1,:], kind='linear')
	t   = np.linspace(data_tr[0,0],data_tr[0,-1], data_tr.shape[1])
	r_t = F_rt(t)

	#plt.plot(t,r_t)
	#plt.show()
	
	C_t = ComputeCt(r_t)
	C_t /= C_t[0]
	plt.plot(t[0:len(C_t)/20] - t[0], C_t[0:len(C_t)/20])
	plt.xlabel('fs')
	plt.ylabel('Correlation')
	plt.show()
	

if __name__=="__main__":
	main()
