import numpy as np
import matplotlib.pyplot as plt
import os

dir = '../result/'

filenames=os.listdir(dir)

B = np.loadtxt(dir+'Heart')

for filename in filenames:
    if filename.startswith('Heart'): continue
    plt.cla()
    A = np.loadtxt(dir+filename)
    if filename.startswith('problemA('):
        plt.plot(A[0,:],A[1,:],label = 'f' )
        plt.plot(A[0,:],A[2,:],label = 'spline(N=6)' )
        plt.plot(A[0,:],A[3,:],label = 'spline(N=11)' )
        plt.plot(A[0,:],A[4,:],label = 'spline(N=21)' )
        plt.plot(A[0,:],A[5,:],label = 'spline(N=41)' )
        plt.plot(A[0,:],A[6,:],label = 'spline(N=81)' )
    if filename.startswith('problemAError'):
        plt.plot(A[0,:],A[1,:],label = 'complete' )
        plt.plot(A[0,:],A[2,:],label = 'nature' )
        plt.plot(A[0,:],A[3,:],label = 'second' )
        plt.plot(A[0,:],A[4,:],label = 'notAKnot' )
        plt.plot(A[0,:],A[5,:],label = 'periodic' )
        if filename.startswith('problemAErrors('):
            for i in range(1,6):
                i = int(i)
                x = A[0,:]
                y = A[i,:]
                for a,b in zip(x,y):
                    plt.text(a, b+0.05, '%.0f' % b, ha='center', va= 'bottom',fontsize=7)
        plt.yscale('log')
        plt.xlabel('N')
        plt.ylabel('error')
    if filename.startswith('problemBCD'):
        plt.plot(A[:,0],A[:,1],label = 'f' )
        plt.plot(A[:,0],A[:,2],label = 'CardinalSplines<2>')
        plt.plot(A[:,0],A[:,3],label = 'CardinalSplines<3>')
    if filename.startswith('problemE'):
        plt.plot(A[:,0],A[:,1],label = 'fitted by spline' )
        plt.plot(B[:,0],B[:,1],label = 'Heart(exact)', color = 'r', linestyle = '--' )
    plt.title(filename)
    plt.legend()
    plt.savefig('../fig/'+filename+'.jpg')