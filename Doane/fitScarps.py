import os
os.environ["OMP_NUM_THREADS"] = '1'
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pdb
from sklearn import mixture
from scipy.optimize import minimize
from scipy.linalg import inv
#import warnings
#warnings.filterwarnings('ignore')


def fitInitCondit(y0):
    zInitCondit = np.zeros_like(x)
    scarpSlope = y0 - ruptureSlope*x
    zInitCondit = np.copy(scarpSlope)
    zInitCondit[scarpSlope>footFarf]=footFarf[scarpSlope>footFarf]
    zInitCondit[scarpSlope<hangFarf]=hangFarf[scarpSlope<hangFarf]
    mass = np.trapz(zInitCondit - y)**2
    return(mass)

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

def evolScarp(ZTest, p):
    ## Evolve scarp for testing ##
    E0 = 0.0
    E1, L0 = p

    t=0
    Z0 = np.copy(ZTest)
    while t<T:
        slp = np.gradient(Z0)/dx
        mat= convMatMake(Z0, E1, L0)
        Q = np.sum(mat,axis=0)
        dz = -np.gradient(Q)/dx
        Z0 += dz
        t += dt
        if dt>T:
            print(dt)
    return(Z0)

def GaussIterate(p1, p2, s1, s2, s3, damper):
    B = np.array([p1[0], p2[0]])
    Jacob = np.array(([(s2-s1)/(p1[1] - p1[0])], [(s3 - s1)/(p2[1]-p2[0])]))
    B = np.mat(B).T
    J = np.mat(Jacob).T
    S = np.mat(s1).T
    damping = damper*np.mat(np.array(([1, 0], [0, 1])))

    Bnew = B - inv(J.T.dot(J)+ damping).dot(J.T.dot(S))
    return(Bnew)

#print(os.path.expanduser('~'))
pName = os.path.expanduser('~') + '\DOI\Gray, Harrison J - NonLocalScarps\GIS Scarp Tracing'
fName = '\scarp_datasets_processed.xlsx'

dataset = pd.read_excel(pName + fName)
profileNames = dataset.scarp_unique_ID.unique()

goodScarps = np.load('goodProfiles.npy', allow_pickle = True)

collectData = pd.read_csv(pName + '\\nonlocalScarpParams.csv')
pdb.set_trace()

## Gaussian mixture model for isolating parts of the scarp
G1 = mixture.GaussianMixture(n_components=3, covariance_type='full')
#G4 = mixture.GaussianMixture(n_components=3, covariance_type='spherical')

#collectData = pd.DataFrame(columns= ['name', 'age', 'E1', 'L0', 'Sum of Errors', 'Converged'])
noFit = {'name': []}
noFit = pd.DataFrame(noFit)

goodProfiles = []
ruptureSlope = np.tan(70*np.pi/180)
n = 1
N = len(goodScarps)
for name in goodScarps:
    if name in collectData['name'].values:
        print('skipping, already done')
        continue

    temp = dataset.loc[dataset['scarp_unique_ID']==name]
    x = temp['x_spline'].to_numpy()
    x-=np.mean(x)
    dx = x[1]-x[0]
    y = temp['y_spline'].to_numpy()
    y-=np.mean(y)
    age = temp['age'].to_numpy()[0]*1000

    if age>50000:
        print('skipping, age too old')
        continue

    slp = np.gradient(y,x)
    slope = np.abs(np.arctan(slp)*180/np.pi)
    obs = np.array(np.hstack((np.mat(x).T, np.mat(np.abs(slp)).T)))
    #obs = np.array(np.mat(slp).T)

    G1.fit(obs)
    labs = G1.predict(obs)
    
    print(name, n/N, age)
    n+=1
    # Get unique identifiers for each part of the scarp
    tempLabs = [0, 1, 2]
    fLab = labs[0]
    tempLabs.remove(labs[0])
    hLab = labs[-1]
    tempLabs.remove(labs[-1])
    sLab = tempLabs[0]

    zFWall = y[labs==fLab]
    xFWall = x[labs==fLab]

    zHWall = y[labs==hLab]
    xHWall = x[labs==hLab]

    zScarp = y[labs == sLab]
    xScarp = x[labs ==sLab]

    slpFWall = slp[labs==fLab]
    slpHWall = slp[labs==hLab]
    slpScarp = slp[labs==sLab]

    # #slopes = G1.medians_[:,1]
    centers = G1.means_[:,0]

    fSlope = np.median(slpFWall)
    fCenter = centers[fLab]
    hSlope = np.median(slpHWall)
    hCenter = centers[hLab]

    hIntercept = np.mean(zHWall)- hSlope*np.mean(xHWall)
    fIntercept = np.mean(zFWall) - fSlope*np.mean(xFWall)

    footFarf = fSlope*x + fIntercept
    hangFarf = hSlope*x + hIntercept

    yInt = minimize(fitInitCondit, 5)
    zInitCondit = np.zeros_like(x)
    try:
        scarpSlope = yInt.x[0] - ruptureSlope*x
    except:
        continue
    zInitCondit = np.copy(scarpSlope)
    zInitCondit[scarpSlope>footFarf]=footFarf[scarpSlope>footFarf]
    zInitCondit[scarpSlope<hangFarf]=hangFarf[scarpSlope<hangFarf]


    T = age
    Sc = 1.5

    t = 0 
    dt = 0.5
    switch = 1 #Switch for when tolerance drops below the desired value
    damper = 0.5
    sBest = 100.
    sOld = 100.

    k = 0

    dE = 0.0001
    dL = 0.001

    ## Initial guesses for Entrainment rate and Lambda ##
    #EP = [0.001] #Entrainment
    #EL = [0.1] #Lambda
    B = np.mat([0.001, 0.1]).T #Column Vector of initial guesses
        
    while switch==1:
        EP = [B[0,0], B[0,0]+dE] 
        EL = [B[1,0], B[1,0]+dL]

        ## Evolve scarp according to three combinations of parameters to build the Jacobian matrix
        z1 = evolScarp(zInitCondit, [EP[0], EL[0]]) #First Guesses
        s1 = y[labs==sLab] - z1[labs==sLab] #Evaluate residuals 
        z2 = evolScarp(zInitCondit, [EP[1], EL[0]]) #Second Guess and First Guess
        s2 = y[labs==sLab] - z2[labs==sLab] 
        z3 = evolScarp(zInitCondit, [EP[0], EL[1]]) #First Guess and Second Guess
        s3 = y[labs==sLab] - z3[labs==sLab]

        BNew = GaussIterate(EP, EL, s1, s2, s3, damper) #Matrix multiplication steps to determine array of new guesses
        dS = sOld - np.sum(s1**2)

        if dS>=0:
            
            damper/=2
        else:
            damper*=2
        
        if np.sum(s1**2)<sOld:
            BBest = np.copy(B)
            sBest = np.copy(np.sum(s1**2))

        sOld = np.sum(s1**2)
      
        if any(BNew<0) or any(BNew>1):
            dB = BNew - B #Change in guesses
            print('negative values')
            while any(BNew<0) or any(BNew>2):
                dB /=2
                BNew = B + dB
        dB = BNew - B
        print(k, B, np.sqrt(dB[0]**2 + dB[1]**2), sOld) 

        k+=1
        if k>=15:
            noFit.loc[len(noFit)]={'name': name}
            temp = {'name': name, 'age': age, 'E1': BBest[0,0], 'L0': BBest[1,0], 'Sum of Errors': sBest, 'Converged': 'No'}
            collectData.loc[len(collectData)] = temp
            collectData.to_csv(pName + '/nonlocalScarpParams.csv', index= False)
            switch = 0
            continue

        #Evaluate if change was small enough to exit loop and call it good!
        if (np.sqrt(dB[0]**2 + dB[1]**2)<0.005 and sOld==sBest):
            switch = 0
            temp = {'name': name, 'age': age, 'E1': B[0,0], 'L0': B[1,0], 'Sum of Errors': sOld, 'Converged': 'Yes'}
            #temp = pd.DataFrame([[name, age, BNew[0, 0], BNew[1, 0], sOld]], columns= ['name', 'age', 'E1', 'L0', 'Sum of Errors'])
            collectData.loc[len(collectData)] = temp
            collectData.to_csv(pName + '/nonlocalScarpParams.csv', index= False)

        B = BNew
        if sOld<sBest:
            sBest = sOld
    plt.scatter(x, y, c = labs)
    plt.plot(x, zInitCondit, '-k')
    plt.plot(x, z1, '-k')
    plt.pause(5)
    plt.close()

collectData.to_csv('scarpParameters.csv', sep = ',')
noFit.to_csv('badFitScarps.csv', sep = ',')

pdb.set_trace()

