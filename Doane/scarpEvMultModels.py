import os
os.environ["OMP_NUM_THREADS"] = '1'
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pdb
from sklearn import mixture
from scipy.optimize import minimize
from scipy.linalg import inv

def convMatMake(Z0, E1, L0, Sc):
    #disent = 2*(Sc - np.abs(slp))
    slp = np.gradient(Z0)/dx
    disent = (Sc + slp)/(Sc - slp)*dx/L0
    disent[disent<0]=0.

    mat = np.zeros((len(slp), len(slp)))

    p_x = (1./2.)*(1.0 - slp/Sp)
    p_x[-slp>Sp] = 1.0
    p_x[-slp<-Sp] = 0.0
    n_x = 1.0 - p_x

    p_x[:] = 1
    for i in range(len(slp)):
        #ESlp = np.random.exponential((E1*np.abs(slp[i])))
        ESlp = E1*np.abs(slp[i])
        E = (ESlp)*dt*dx
        dis1 = np.cumsum(disent[i:])
        dis2 = dis1[-1] + np.cumsum(disent[:i])
        mat[i,i:]=p_x[i]*E*np.exp(-dis1)
        mat[i,:i] =p_x[i]*E*np.exp(-dis2)
    #pdb.set_trace()
    return(mat)

dx = 0.1
L1 = 0.05
E1 = 0.001
E2 = 0.001
L2 = 0.5

Sc = 1.5
Sp = 0.001
dt = 1

farSlope = -0.00
scarpSlope = -2.0

x = np.arange(-10, 10, dx)
fWall = x*farSlope + 3.
hWall = x*farSlope - 3.
z = x*scarpSlope

z[z>fWall] = fWall[z>fWall]
z[z<hWall] = hWall[z<hWall]

#z[x>3]+=.2

#pdb.set_trace()

T = 1000
t = 0
alpha = 0.01
while t<T:
    matOne = convMatMake(z, E1, L1, 1.0)
    matTwo = convMatMake(z, E2, L2, 1.0)
    #matNeg = convMatMake(np.flip(z), E1, L0)
    QOne = np.sum(matOne, axis=0)
    QTwo = np.sum(matTwo, axis=0)
    Q = QOne + QTwo
    dz = -np.gradient(Q)/dx

    #pdb.set_trace()
    z+=dz
    t+=dt
    if np.remainder(t,100)==0:
        plt.plot(x[np.abs(x)<10],z[np.abs(x)<10],'-k')

plt.axis('equal')
plt.show(block = False)
pdb.set_trace()


