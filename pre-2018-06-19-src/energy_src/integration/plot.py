import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

data = np.loadtxt("discrete.txt")

plt.plot(data[:,0],data[:,1],'ro')

data = np.loadtxt("spline.txt")

plt.plot(data[:,0],data[:,1],'k.')

print(simps(data[5:,1],data[5:,0]))

plt.show()
