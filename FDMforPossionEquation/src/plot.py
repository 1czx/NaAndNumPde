import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os 

def plot_heatmap(solution, save_path):
    plt.imshow(solution, cmap='hot', origin='lower')
    plt.colorbar()
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(save_path)
    plt.savefig(outputDir+save_path+'.jpg')
    plt.close()

dir = '../result/'

filenames=os.listdir(dir)

outputDir = '../fig/'

for filename in filenames:
    if filename.endswith('E'):
        solution = np.loadtxt(dir+filename)
        plot_heatmap( solution, filename )

x = [1.0/8,1.0/16,1.0/32,1.0/64]

x = np.array(x)

x = np.log2(x)

for filename in filenames:
    if filename.endswith('Norm'):
        A = np.loadtxt(dir+filename)
        A = np.log2(A)
        plt.plot(x,A[:,0],label="L1-Norm")
        plt.plot(x,A[:,1],label="L2-Norm")
        plt.plot(x,A[:,2],label="LInf-Norm")
        plt.xlabel('log(h)')
        plt.ylabel('log(||E||)')
        plt.title(filename)
        plt.legend()
        plt.savefig(outputDir+filename+'jpg')
        plt.cla()

def test1(x,y):
    return np.power(np.e,y+np.sin(x))
def test2(x,y):
    return np.sin(np.pi*x)*np.sin(np.pi*y)
def test3(x,y):
    return x**3+y**3

x = np.linspace(0, 1, 100)
y = np.linspace(0, 1, 100)
X, Y = np.meshgrid(x, y)
Z = test1(X,Y)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z (exp(y+sin(x)))')
plt.title('Test1True')
plt.savefig(outputDir+'Test1True.jpg')
plt.close()

Z = test2(X,Y)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z (sin(pi*x)*sin(pi*y))')
plt.title('Test2True')
plt.savefig(outputDir+'Test2True.jpg')
plt.close()

Z = test3(X,Y)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z (x^3+y^3)')
plt.title('Test3True')
plt.savefig(outputDir+'Test3True.jpg')
plt.close()