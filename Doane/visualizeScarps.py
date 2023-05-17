import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pdb
from sklearn import mixture
import os

#print(os.path.expanduser('~'))
pName = os.path.expanduser('~') + '\DOI\Gray, Harrison J - NonLocalScarps\GIS Scarp Tracing'
fName = '\scarp_datasets_processed.xlsx'

dataset = pd.read_excel(pName + fName)
profileNames = dataset.scarp_unique_ID.unique()

## Gaussian mixture model for isolating parts of the scarp
G1 = mixture.GaussianMixture(n_components=2, covariance_type='full')
G2 = mixture.GaussianMixture(n_components=3, covariance_type='tied')
G3 = mixture.GaussianMixture(n_components=3, covariance_type='diag')
G4 = mixture.GaussianMixture(n_components=2, covariance_type='spherical')

goodProfiles = []
for name in profileNames:
    temp = dataset.loc[dataset['scarp_unique_ID']==name]
    x = temp['x_spline']
    y = temp['y_spline']

    slp = np.gradient(y,x)
    slope = np.arctan(slp)
    #obs = np.array(np.hstack((np.mat(x).T, np.mat(np.abs(slp)).T)))
    obs = np.array(np.mat(slp).T)

    G4.fit(obs)
    labs = G4.predict(obs)

    plt.scatter(x, y, c = labs)
    plt.show(block = False)
    
    keep = input('Does the scarp look good? y/n')
    if keep=='y':
        goodProfiles = np.append(goodProfiles, name)

    plt.close()

pdb.set_trace()
np.save('goodProfiles.npy', goodProfiles, allow_pickle=True)
