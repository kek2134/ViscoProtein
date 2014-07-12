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

def ComputeJw(x_t,dt):
	dx_t = x_t - np.mean(x_t)
	x_w = fft.fft(dx_t)
	J_w = np.absolute(x_w * np.conjugate(x_w))
        w = fft.fftfreq(len(x_t), dt)
        J_w *= w
        print "Max of J_w: {}".format(max(J_w))
        point_density = len(w) / (max(w) - min(w))
        print "Point density: {}, Number of points: {}, Frequency Axis: {}".format(point_density, len(w), max(w)-min(w))
        J_w /= point_density

        return J_w, w


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
	
        run_Ct = True
        if run_Ct:
	    C_t = ComputeCt(r_t)
	    C_t /= C_t[0]
            fraction = 5
	    plt.plot(t[0:len(C_t)/fraction] - t[0], C_t[0:len(C_t)/fraction])
	    plt.xlabel('fs')
	    plt.ylabel('Correlation')
	    plt.show()

	    #plt.loglog(t[0:len(C_t)/fraction] - t[0], -np.log(C_t[0:len(C_t)/fraction]))
	    #plt.xlabel('fs')
	    #plt.ylabel('Correlation')
	    #plt.show()


        run_Jw = False
        if run_Jw:
            
            J_w, w = ComputeJw(r_t, t[1]-t[0])
    
            J_w = J_w[:len(w)/2]
            w   =   w[:len(w)/2]
    
            D = 500
            J_w_red = np.zeros(len(J_w)/D)
            w_red = np.zeros(len(J_w)/D)
            for i in xrange(len(J_w)/D):
                J_w_red[i] = np.mean(J_w[i*D:(i*D)+D])
                w_red[i] = w[i*D]
    
    
            plt.loglog(w_red,J_w_red)
            plt.xlabel('1/ps or something')
            plt.ylabel('r^2 or something')
            plt.show()
    
	

if __name__=="__main__":
	main()
