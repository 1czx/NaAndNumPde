import numpy as np
import matplotlib.pyplot as plt
import os 

dir = '../result/'

filenames=os.listdir(dir)

outputDir = '../fig/'

for filename in filenames:
    A = np.loadtxt(dir+filename)
    if filename.endswith('A'):
        fig, ax1 = plt.subplots()
        color = 'tab:blue'
        ax1.set_xlabel('log(num of steps)')
        ax1.set_ylabel('log(||E||)')
        ax1.plot(np.log2(A[:,0]), np.log2(A[:,1]),color = color)
        ax1.tick_params(axis='y', labelcolor = color)
        color = 'tab:red'
        ax2 = ax1.twinx()
        ax2.set_ylabel('log(Time/ms)')
        ax2.plot(np.log2(A[:,0]), np.log(A[:,2]),color = color)
        ax2.tick_params(axis='y',labelcolor = color)
        # fig.tight_layout()
        # plt.plot(np.log2(A[:,0]),np.log2(A[:,1]))
        # plt.xlabel('log(num of steps)')
        # plt.ylabel('log(||E||)')
        plt.title(filename)
        plt.savefig(outputDir+filename+'.jpg')
        plt.cla()
        plt.close()
    if filename.endswith('P'):
        plt.plot(A[:,0],A[:,1])
        plt.title(filename)
        plt.savefig(outputDir+filename+'.jpg')
        plt.cla()