import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os 

def read_line_to_array(file_path, line_number):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        if line_number <= len(lines):
            line = lines[line_number - 1]
            data = np.array([float(value) for value in line.split()])
            return data

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

x = [1.0/32,1.0/64,1.0/128,1.0/256]
n = [32,64,128,256]

x = np.array(x)
n = np.array(n)

x = np.log2(x)
n=np.log2(n)

for filename in filenames:
    if filename.endswith('Norm'):
        A = np.loadtxt(dir+filename)
        A = np.log2(A)
        plt.plot(x,A,label="LInf-Norm")
        plt.xlabel('log(h)')
        plt.ylabel('log(||E||)')
        plt.title(filename)
        plt.legend()
        plt.savefig(outputDir+filename+'.jpg')
        plt.cla()
    if filename.endswith('E'):
        if filename.startswith('2'):
            solution = np.loadtxt(dir+filename)
            plot_heatmap( solution, filename )
        else:
            A = np.loadtxt(dir+filename,)
            plt.plot(A[:,0],A[:,1])
            plt.xlabel('x')
            plt.ylabel('AbsError')
            plt.title(filename)
            plt.savefig(outputDir+filename+'.jpg')
            plt.cla()
    if filename.endswith('Res'):
        plt.ylim(0,1)
        A = read_line_to_array(dir+filename,1)
        k = np.array(range(2,A.size))
        A = A[2:A.size]/A[1:A.size-1]
        plt.plot(k,A,label='n=32')
        A = read_line_to_array(dir+filename,2)
        k = np.array(range(2,A.size))
        A = A[2:A.size]/A[1:A.size-1]
        plt.plot(k,A,label='n=64')
        A = read_line_to_array(dir+filename,3)
        k = np.array(range(2,A.size))
        A = A[2:A.size]/A[1:A.size-1]
        plt.plot(k,A,label='n=128')
        A = read_line_to_array(dir+filename,4)
        k = np.array(range(2,A.size))
        A = A[2:A.size]/A[1:A.size-1]
        plt.plot(k,A,label='n=256')
        plt.xlabel('iter times')
        plt.ylabel('||r_k+1||/||r_k||)')
        plt.title(filename)
        plt.legend()
        plt.savefig(outputDir+filename+'.jpg')
        plt.cla()
    if filename.endswith('Time'):
        A = np.loadtxt(dir+filename)
        A = np.log2(A)
        plt.plot(n,A)
        plt.xlabel('log(n)')
        plt.ylabel('log(||Time||)')
        plt.title(filename)
        plt.savefig(outputDir+filename+'.jpg')
        plt.cla()
    
    
def test12D(x,y):
    return np.power(np.e,y+np.sin(x))
def test22D(x,y):
    return np.sin(np.pi*x)*np.sin(np.pi*y)
def test32D(x,y):
    return x**3+y**3

def test11D(x):
    return np.power(np.e,x+np.sin(x))
def test21D(x):
    return np.sin(np.pi*x)
def test31D(x):
    return x**3

x = np.linspace(0, 1, 100)

y = test11D(x)
plt.plot(x,y)
plt.xlabel('X')
plt.ylabel('Y(exp(x+sin(x)))')
plt.title("1DTest1Ture")
plt.savefig(outputDir+'1DTest1Ture.jpg')
plt.cla()

y = test21D(x)
plt.plot(x,y)
plt.xlabel('X')
plt.ylabel('Y(sin(pi*x))')
plt.title("1DTest2Ture")
plt.savefig(outputDir+'1DTest2Ture.jpg')
plt.cla()

y = test31D(x)
plt.plot(x,y)
plt.xlabel('X')
plt.ylabel('Y(x^3)')
plt.title("1DTest3Ture")
plt.savefig(outputDir+'1DTest3Ture.jpg')
plt.cla()

y = np.linspace(0, 1, 100)
X, Y = np.meshgrid(x, y)
Z = test12D(X,Y)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z (exp(y+sin(x)))')
plt.title('2DTest1True')
plt.savefig(outputDir+'2DTest1True.jpg')
plt.close()

Z = test22D(X,Y)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z (sin(pi*x)*sin(pi*y))')
plt.title('2DTest2True')
plt.savefig(outputDir+'2DTest2True.jpg')
plt.close()

Z = test32D(X,Y)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z (x^3+y^3)')
plt.title('2DTest3True')
plt.savefig(outputDir+'2DTest3True.jpg')
plt.close()