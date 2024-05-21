import numpy as np
import matplotlib.pyplot as plt

A = np.loadtxt('test.txt')

plt.plot(A[:,0],A[:,1])
plt.savefig('1.jpg')