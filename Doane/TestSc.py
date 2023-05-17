import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import mixture
from scipy.optimize import minimize
from scipy.sparse import eye
from numpy.linalg import inv
import pdb

def convMatMake(Z0, E1, L0):
    #disent = 2*(Sc - np.abs(slp))
    slp = np.gradient(Z0)/dx
    disent = (Sc - np.abs(slp))/(Sc + np.abs(slp))*dx/L0
    disent[disent<0]=0.
    disent[slp>0]=1.

    mat = np.zeros((len(slp), len(slp)))
    indTemp = np.arange(len(slp))
    #dt = dx/np.sqrt(np.max(np.abs(slp)))
    #if dt>1:
    #    dt=1.
    dt = 1/2
    for i in range(len(slp)):
        #ESlp = np.random.exponential((E1*np.abs(slp[i])))
        ESlp = E1*(np.abs(slp[i]))**(2/2)
        #if np.max(ESlp)>0.01:
        #    dt = 0.01/np.max(ESlp)
        #else:
        #    dt=1.
        #ESlp[slp>0]=0
        E = (ESlp)*dt*dx
        dis1 = np.cumsum(disent[i:])
        dis2 = dis1[-1] + np.cumsum(disent[:i])
        mat[i,i:]=E*np.exp(-dis1)
        mat[i,:i] =E*np.exp(-dis2)
    #pdb.set_trace()
    return(mat, dt)

def evolScarp(ZTest, p):
    ## Evolve scarp for testing ##
    E0 = 0.0
    E1, L0 = p

    t=0
    #T=500
    Z0 = np.copy(ZTest)
    while t<T:
        slp = np.gradient(Z0)/dx
        mat, dt = convMatMake(Z0, E1, L0)
        Q = np.sum(mat,axis=0)
        dz = -np.gradient(Q)/dx
        Z0 += dz
        t += dt
        if dt>T:
            print(dt)
    return(Z0)

##Transport and Scarp Parameters##

L = 2.
slope = np.tan(80*np.pi/180)
H = slope*L

dx = 0.1
x = np.arange(-10*L, 10*L, dx)
zInit = np.zeros_like(x)
zInit = -slope*x
zInit[x<-L]=H - 0.2*x[x<-L]
zInit[x>L]=-H - 0.1*x[x>L]

T = 250

E1 = 0.001
L0 = 0.2
p = [E1, L0]

talusSlp = []
testSlp = np.array([0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
for Sc in testSlp:
    T=100
    print(Sc)
    z = evolScarp(zInit, p)
    plt.plot(x[np.abs(x)<4], z[np.abs(x)<4])
    pts = plt.ginput(2,timeout=-1)
    slpTemp = (pts[0][1] - pts[1][1])/(pts[0][0]-pts[1][0])
    talusSlp = np.append(talusSlp, np.abs(slpTemp))
    plt.close()

plt.plot(testSlp ,talusSlp/testSlp, 'ok')
plt.show(block = False)
pdb.set_trace()