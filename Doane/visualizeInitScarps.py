import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pdb
from sklearn import mixture
from scipy.optimize import minimize
import os

def fitInitCondit(y0):
    zInitCondit = np.zeros_like(x)
    scarpSlope = y0 - ruptureSlope*x
    zInitCondit = np.copy(scarpSlope)
    zInitCondit[scarpSlope>footFarf]=footFarf[scarpSlope>footFarf]
    zInitCondit[scarpSlope<hangFarf]=hangFarf[scarpSlope<hangFarf]
    mass = np.trapz(zInitCondit - y)**2
    return(mass)


#print(os.path.expanduser('~'))
pName = os.path.expanduser('~') + '\DOI\Gray, Harrison J - NonLocalScarps\GIS Scarp Tracing'
fName = '\scarp_datasets_processed.xlsx'

dataset = pd.read_excel(pName + fName)
profileNames = dataset.scarp_unique_ID.unique()

goodScarps = np.load('goodProfiles.npy', allow_pickle = True)

## Gaussian mixture model for isolating parts of the scarp
G1 = mixture.GaussianMixture(n_components=3, covariance_type='full')
#G4 = mixture.GaussianMixture(n_components=3, covariance_type='spherical')

goodProfiles = []
ruptureSlope = np.tan(70*np.pi/180)
for name in profileNames:
    temp = dataset.loc[dataset['scarp_unique_ID']==name]
    x = temp['x_spline'].to_numpy()
    x-=np.mean(x)
    y = temp['y_spline'].to_numpy()
    y-=np.mean(y)

    slp = np.gradient(y,x)
    slope = np.abs(np.arctan(slp)*180/np.pi)
    obs = np.array(np.hstack((np.mat(x).T, np.mat(np.abs(slp)).T)))
    #obs = np.array(np.mat(slp).T)

    G1.fit(obs)
    labs = G1.predict(obs)
    
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

    plt.scatter(x, y, c = labs)
    plt.plot(x, zInitCondit, '-k')
    plt.show(block = False)
    
    keep = input('Does the scarp look good? y/n')
    if keep=='y':
        goodProfiles = np.append(goodProfiles, name)
        np.save('goodProfiles.npy', goodProfiles, allow_pickle=True)

    plt.close()
   

pdb.set_trace()

