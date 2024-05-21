import numpy as np
import matplotlib.pyplot as plt
import os 

dir = '../result/'

filenames=os.listdir(dir)

outputDir = '../fig/'

for filename in filenames:
    A = np.loadtxt(dir+filename)
    if filename.endswith('A'):
        plt.plot(np.log2(A[:,0]),np.log2(A[:,1]))
        plt.title(filename)
        plt.savefig(outputDir+filename+'.jpg')
        plt.cla()
    if filename.endswith('P'):
        plt.plot(A[:,0],A[:,1])
        plt.title(filename)
        plt.savefig(outputDir+filename+'.jpg')
        plt.cla()