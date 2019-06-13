import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("discrete.txt")

plt.plot(data[:,0],data[:,1],'o')

data = np.loadtxt("continuous.txt")

plt.plot(data[:,0],data[:,1],'.')


data = np.loadtxt("data/psivsr.txt")

plt.plot(data[:,0],data[:,1],'.')

plt.show()
