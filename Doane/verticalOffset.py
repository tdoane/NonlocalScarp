import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pdb

#pName = 'C:\\Users\\tdoane\\DOI\\Gray, Harrison J - NonLocalScarps\\GIS Scarp Tracing'
fName = 'Datasets.xlsx'
metaName = 'Metadata.csv'

data = pd.read_excel(fName, sheet_name='Data')
mData = pd.read_excel(fName, sheet_name='Metadata')

source = data['SourcePub'].unique()
goodScarps= []
goodPubs = []
offset = []
ages = []
skew = []
dx = 0.1
for pub in source:
    temp = data.loc[data["SourcePub"]==pub]
    profiles = temp['Profile_name'].unique()
    mTemp = mData.loc[mData["SourcePub"]==pub]
    for name in profiles:
        name = str(name)
        prof = temp.loc[temp['Profile_name']==name]
        aData = mTemp.loc[mTemp['Profile_name'].astype('str')==name]
        age = aData['age'].to_numpy()[0]
        x = prof['X'].to_numpy()
        z = prof['Y'].to_numpy()
        xTemp, inds = np.unique(x, return_index=True)
        x = x[inds]
        z = z[inds]
        print(pub, name, age)
        plt.plot(x,z, 'ok')
        plt.plot(x,z, '-k')
        plt.show(block = False)
        keeper = input('Is the scarp morphologically clear? y/n')
        if keeper=='n':
            print('skipping unusual scarp form')
            plt.close()
            continue    
        
        print('Now click two points along the far-field slope on the footwall (upslope side of the scarp). First click the left-most point, then just upslope of the scarp, then press enter. The first point that you select will appear as a red cross, the second one will not')
        XYfoot = plt.ginput(2)
        print('Now repeat this on the hangingwall (downslope of the scarp)')
        XYhang = plt.ginput(2)
        plt.close()

        footSlope = (XYfoot[1][1] - XYfoot[0][1])/(XYfoot[1][0]- XYfoot[0][0])
        hangSlope = (XYhang[1][1] - XYhang[0][1])/(XYhang[1][0]- XYhang[0][0])
        
        footIntercept = XYfoot[0][1] - XYfoot[0][0]*footSlope
        hangIntercept = XYhang[0][1] - XYhang[0][0]*hangSlope

        xEqual = np.arange(min(x), max(x), dx) 
        zEqual = np.interp(xEqual, x, z)

        ind = np.where((xEqual>=XYfoot[1][0]) & (xEqual<=XYhang[0][0]))
        xScarp = xEqual[ind]
        zScarp = zEqual[ind]
       
        plt.plot(xScarp,zScarp, 'ok')
        plt.plot(x, footIntercept + footSlope*x)
        plt.plot(x, hangIntercept + hangSlope*x)
        plt.grid()
        plt.show(block = False)
        
        print('click two points along the new lines for the vertical offset of the scarp.')
        offsetPoints = plt.ginput(2)
        offsetHeight = np.abs(offsetPoints[0][1] - offsetPoints[1][1])
        plt.close()

        xTemp = xScarp/offsetHeight
        xTemp-= np.mean(xTemp)
        zTemp = zScarp/offsetHeight
        #zTemp = (zTemp+s0*xTemp)
        slpTemp = np.abs(np.gradient(zTemp, xTemp))
        norm = np.trapz(np.abs(slpTemp),xTemp)
        moment1 = np.trapz(xTemp*np.abs(slpTemp), xTemp)
        moment2 = np.trapz((xTemp-moment1)**2*np.abs(slpTemp), xTemp)#/np.trapz(slpTemp,x)
        moment3 = np.trapz(((xTemp-moment1)/np.sqrt(moment2))**3*np.abs(slpTemp), xTemp)#/np.trapz(slpTemp,x)

        print(moment3)
        plt.plot(xTemp, slpTemp, '-k')
        plt.show(block = False)
        check = input('Does this look good? y/n ')
        if check=='n':
            continue  
        else:
            plt.close()

        offset.append(offsetHeight)
        goodScarps.append(name)
        goodPubs.append(pub)
        ages.append(age)
        skew.append(moment3)

        next = input("Do you want to continue? y/n ")
        if next == "n":
            scarpCollect = {'Source': goodPubs, 'Profiles': goodScarps, 'Offset': offset, 'Asymmetry': skew, 'Age': ages}
            np.save('scarpOffsetAgeSkew.npy', scarpCollect, allow_pickle=True)
            pdb.set_trace()
            break