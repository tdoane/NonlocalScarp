import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import pandas as pd
from numpy.linalg import inv
from multiprocessing import Pool
from scipy.optimize import curve_fit, minimize
from sklearn import mixture
import pdb

def convMatMake(Z0, E1, L0):
    #disent = 2*(Sc - np.abs(slp))
    slp = np.gradient(Z0)/dx
    disent = (Sc - np.abs(slp))/(Sc + np.abs(slp))*dx/L0
    disent[disent<0]=0.

    mat = np.zeros((len(slp), len(slp)))
    indTemp = np.arange(len(slp))
    for i in range(len(slp)):
        #ESlp = np.random.exponential((E1*np.abs(slp[i])))
        ESlp = E1*np.abs(slp[i])
        E = (ESlp)*dt*dx
        dis1 = np.cumsum(disent[i:])
        dis2 = dis1[-1] + np.cumsum(disent[:i])
        mat[i,i:]=E*np.exp(-dis1)
        mat[i,:i] =E*np.exp(-dis2)
    #pdb.set_trace()
    return(mat)

def evolScarp(Z, p):
    ## Evolve scarp for testing ##
    E0 = 0.0
    Sc = 1.2
    E1, L0 = p

    t=0
    #T=500
    Z0 = np.copy(Z)
    while t<T:
        slp = np.gradient(Z0)/dx
        mat = convMatMake(Z0, E1, L0)
        Q = np.sum(mat,axis=0)
        dz = -np.gradient(Q)/dx
        Z0 += dz
        t += dt
    return(Z0)

def GaussIterate(p1, p2, s1, s2, s3):
    B = np.array([p1[0], p2[0]])
    Jacob = np.array(([(s2-s1)/(p1[1] - p1[0])], [(s3 - s1)/(p2[1]-p2[0])]))
    B = np.mat(B).T
    J = np.mat(Jacob).T
    S = np.mat(s1).T
    Bnew = B - inv(J.T.dot(J)).dot(J.T.dot(S))
    return(Bnew)

def InitCondit(x, center, height, hangSlp, footSlp):
    L = height/slope
    z0 = np.zeros_like(x)
    z0 = -slope*(x-center)
    z0[(x-center)<-L]=height + x[(x-center)<-L]*footSlp
    z0[(x-center)>L] = -height + x[(x-center)>L]*hangSlp
    #z0+=farFSlp*(x-center)
    return(z0)

def fitInitCondit(y0):
    zInitCondit = np.zeros_like(x)
    scarpSlope = y0 - slope*x
    zInitCondit = scarpSlope
    zInitCondit[scarpSlope>footFarf]=footFarf[scarpSlope>footFarf]
    zInitCondit[scarpSlope<hangFarf]=hangFarf[scarpSlope<hangFarf]
    mass = np.trapz(zInitCondit - zTest)**2
    return(mass)

##Transport and Scarp Parameters##

Sc = 1.2
L = 2.
slope = np.tan(60*np.pi/180)
H = slope*L

dx = 0.1
x = np.arange(-10*L, 10*L, dx)
zInit = np.zeros_like(x)
zInit = -slope*x
zInit[x<-L]=H - 0.2*x[x<-L]
zInit[x>L]=-H - 0.1*x[x>L]
#zInit-= x*0.05

## Gaussian mixture model for isolating parts of the scarp
G = mixture.GaussianMixture(n_components=3)

## Evolve scarp for testing ##
dt = 1.
t = 0
T=500

zTest = evolScarp(zInit, [0.005, 0.75])
#zTest += np.random.normal(0, 0.05, size = len(zTest))
slopeTest = np.gradient(zTest)/dx
obs = np.array(np.hstack((np.mat(x).T, np.mat(slopeTest).T)))
G.fit(obs)
labs = G.predict(obs)

# Get unique identifiers for each part of the scarp
fLab = labs[0]
hLab = labs[-1]
sLab = labs[x==0]

zFWall = zTest[labs==fLab]
xFWall = x[labs==fLab]

zHWall = zTest[labs==hLab]
xHWall = x[labs==hLab]

zScarp = zTest[labs == sLab]
xScarp = x[labs ==sLab]

slopes = G.means_[:,1]
centers = G.means_[:,0]

fSlope = slopes[fLab]
fCenter = centers[fLab]
hSlope = slopes[hLab]
hCenter = centers[hLab]

hIntercept = zHWall[0]- hSlope*xHWall[0]
fIntercept = zFWall[-1] - fSlope*xFWall[-1]

footFarf = fSlope*x + fIntercept
hangFarf = hSlope*x + hIntercept

## Fit an Initial condition to the obesrved form that conserves mass ##
yInt = minimize(fitInitCondit, 5)
zInitCondit = np.zeros_like(x)
scarpSlope = yInt.x[0] - slope*x
zInitCondit = scarpSlope
zInitCondit[scarpSlope>footFarf]=footFarf[scarpSlope>footFarf]
zInitCondit[scarpSlope<hangFarf]=hangFarf[scarpSlope<hangFarf]

## 

plt.scatter(x, zTest, c=labs)
plt.plot(x,zInit, '-k')
plt.plot(x, hangFarf, '-r')
plt.plot(x, footFarf, '-y')
plt.plot(x, zInitCondit, '-b')

plt.show(block = False)

## Apply Gauss-Newton Descent Algorithm ##

# Changes in entrainment and L0
dE = 0.0001
dL = 0.01

## Initial guesses for Entrainment rate and Lambda ##
EP = [0.002, 0.002+dE] #Entrainment
EL = [0.5, 0.5+dL] #Lambda
B = np.mat([EP[0], EL[0]]).T #Column Vector of initial guesses

t = 0 
switch = 1 #Switch for when tolerance drops below the desired value
while switch==1:
    EP = [B[0,0], B[0,0]+dE] 
    EL = [B[1,0], B[1,0]+dL]

    ## Evolve scarp according to three combinations of parameters to build the Jacobian matrix
    z1 = evolScarp(zInitCondit, [EP[0], EL[0]]) #First Guesses
    s1 = zTest - z1 #Evaluate residuals 
    z2 = evolScarp(zInitCondit, [EP[1], EL[0]]) #Second Guess and First Guess
    s2 = zTest - z2 
    z3 = evolScarp(zInitCondit, [EP[0], EL[1]]) #First Guess and Second Guess
    s3 = zTest - z3

    BNew = GaussIterate(EP, EL, s1, s2, s3) #Matrix multiplication steps to determine array of new guesses
    dB = BNew - B #Change in guesses
    B = BNew #New Guesses
    print(dB) 
    
    #Evaluate if change was small enough to exit loop and call it good!
    if np.sqrt(dB[0]**2 + dB[1]**2)<0.0001:
        switch = 0


plt.figure()
plt.plot(x, zInitCondit, '-y')
plt.plot(x,z1, '-k')
plt.plot(x, zTest, '-r')
pdb.set_trace()