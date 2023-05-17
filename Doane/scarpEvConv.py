import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.special import erf
import pdb

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))
def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]
def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.

        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp

myColors = np.load('myColors.npy', allow_pickle = 'True').item()
hexlist = myColors['Arizona Sunset']
clrs = plt.get_cmap(get_continuous_cmap(hexlist))

def convMatMake(beta):
    disent =(1 - (np.abs(slp))/Sc)
    disent[disent>1]=1.
    disent[disent<0]=0.
    disent**=beta
    surv = 1 - disent

    mat = np.zeros((len(slp), len(slp)))
    indTemp = np.arange(len(slp))
    for i in range(len(slp)):
        dis1 = np.cumsum(disent[i:]/dx)
        dis2 = dis1[-1] + np.cumsum(disent[:i]/dx)
        #mat[i,i:]= E[i]*np.cumprod(surv[i:])
        mat[i,i:]=E[i]*np.exp(-dis1)
        mat[i,:i] =E[i]*np.exp(-dis2)
    #pdb.set_trace()
    return(mat)


def expMatMake():

    lam = L0*(1./(1 - (np.abs(slp)/Sc)**2))
    lam[np.abs(slp)>Sc]=2*L
    mat = np.zeros((len(slp), len(slp)))
    indTemp = np.arange(len(slp))
    for i in range(len(slp)):
        temp = indTemp-i
        mat[i,:]=E[i]*np.exp(-xArray[temp]/lam[i])
    return(mat)


    #pdb.set_trace()
L = 30.
dx = 0.25
x = np.arange(-L,L,dx)
xArray = np.arange(0, 2*L, dx)
alpha = 3./2.
K = 0.05
T = 10000
t = 0
Sc = 1.2

L0=0.5

E0=0.00001
E1 = 0.0001

k=0
beta = 2

s0 = 0.02

fig, ax = plt.subplots(3,3)

while k<3:
    z=np.ones(len(x))
    z=-3*x
    z[x<-2]=6.
    z[x>2]=-6.
    #pdb.set_trace()
    #z= -2*erf(x)/2.
    z -= s0*x

    t = 0
    T = 10000
    dt=1.

    ax[0, k].plot(x[np.abs(x)<10],z[np.abs(x)<10], '-', color = clrs(0))
    while t<T:
        clr = clrs(float(t/T))
        slp = np.gradient(z)/dx
        print(np.min(slp))
        E = (E0 + E1*np.abs(slp)**(beta))*dt*dx
        mat = expMatMake()
        Q = np.sum(mat,axis=0)
        dz = -np.gradient(Q)/dx
        slpTemp = slp+s0
        if np.remainder(t,100)==0:
            moment1 = np.trapz(x*np.abs(slpTemp), x)
            moment2 = np.trapz(x**2*np.abs(slpTemp), x)#/np.trapz(slpTemp,x)
            moment3 = np.trapz(((x-moment1)/np.sqrt(moment2))**3*np.abs(slpTemp), x)#/np.trapz(slpTemp,x)
            #ax[2,k].plot(t,moment1, 'o', color = clr)
            #ax[3,k].plot(t,moment2, 'o', color = clr)
            ax[2,k].plot(t,moment3, 'o', color = clr)
            if np.remainder(t,500)==0:
                ax[1,k].plot(x[np.abs(x)<10],np.abs(slpTemp[np.abs(x)<10]), '-', color = clr)

        z+=dz
        t+=dt
        print(t)
    beta*=2
    ax[0,k].plot(x[np.abs(x)<10],z[np.abs(x)<10], '-', color = clr)
    ax[0,0].set_aspect('equal')
    ax[0,1].set_aspect('equal')
    ax[0,2].set_aspect('equal')
    k+=1
    #pdb.set_trace()


#moment1 = np.trapz(x*slpTemp, x)#/np.trapz(slpTemp,x)
#moment2 = np.trapz(x**2*slpTemp, x)#/np.trapz(slpTemp,x)
#moment3 = np.trapz(x**3*slpTemp, x)#/np.trapz(slpTemp,x)
#print(moment1, moment2, moment3)

#ax[0].plot(x,z)
#ax[1].plot(x,slpTemp)
plt.show(block = False)
pdb.set_trace()



##Fractional Calculus Model ##
n = 0
N = 1
fig, axs = plt.subplots(2,1)
while n<N:
    z = np.zeros_like(x)
    k = np.fft.fftfreq(len(z))

    count=0
    num = 100.
    T = 1000
    while count<num:
        offset = np.random.weibull(2.,1)*.5
        z[np.abs(x)<100.]+=offset
        ftZ = np.fft.fft(z)
        ftZT=ftZ*np.exp(K*T*(1j*k)**alpha)
        z = np.fft.ifft(ftZT)
        z = np.real(z)
        T = np.random.weibull(3,1)*10000 + 1000.
        #T[T<0] = 1000.
        print(T)

        #axs[0].plot(x[int((L-L/2)/dx):int(2*L/dx)], z[int((L-L/2)/dx):int(2*L/dx)])
        slp = np.gradient(z)/dx

        #axs[1].plot(x[int((L-L/2)/dx):int(2*L/dx)], slp[int((L-L/2)/dx):int(2*L/dx)])
        #axs[0].plot(x,z)
        #axs[1].plot(x,slp)

        count+=1
    axs[0].plot(x,z)
    axs[1].plot(x,slp)
    n+=1
plt.show()
