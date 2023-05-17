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

#myColors = np.load('myColors.npy', allow_pickle = 'True').item()
#hexlist = myColors['Arizona Sunset']
#clrs = plt.get_cmap(get_continuous_cmap(hexlist))
clrs = plt.get_cmap('viridis')

def convMatMake():
    #disent = 2*(Sc - np.abs(slp))
    disent = (Sc - np.abs(slp))/(Sc + np.abs(slp))*dx/L0
    disent[disent<0]=0.

    mat = np.zeros((len(slp), len(slp)))
    indTemp = np.arange(len(slp))
    for i in range(len(slp)):
        #ESlp = np.random.exponential((E1*np.abs(slp[i])))
        ESlp = E1*np.abs(slp[i])
        E = (E0 + ESlp)*dt*dx
        dis1 = np.cumsum(disent[i:])
        dis2 = dis1[-1] + np.cumsum(disent[:i])
        mat[i,i:]=E*np.exp(-dis1)
        mat[i,:i] =E*np.exp(-dis2)
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

L = 20.
dx = 0.1
x = np.arange(-L,L,dx)
xArray = np.arange(0, 2*L, dx)
alpha = 3./2.
K = 0.05
T = 1001
t = 0
Sc = 1.2

E0 = 0.001
E1 = 0.005
L0 = .75

k=0
beta = 1.

s0 = 0.1

fig, ax = plt.subplots(1,3)
slope = 3.
factor = 30
distance = 2.

fault = -x*slope

while k<5:
    #distance+= 0.25
    height = slope*distance
    z=np.ones(len(x))
    z=-slope*x
    z[x<-distance]=height
    z[x>distance]=-height
    #pdb.set_trace()
    #z= -2*erf(x)/2.
    z -= s0*x

    #beta = (1+k)/2.
    beta = 1

    t = 1
    dt=1.
    clr = clrs(float(k/5))

    #ax[0, k].plot(x[np.abs(x)<4*distance],z[np.abs(x)<4*distance], '-', color = clrs(0))
    while t<T:
        eq = np.random.binomial(1, 0.002) ##probabiliyt of a rupture for any time
        slp = np.gradient(z)/dx #calculate the slope of the land surface
        print(np.min(slp)) #Just checking in to make sure it's all reasonable and no numerical instabilities
        if eq==1: #We got a rupture
            offset = 5. #this can be probabilistic in a future model, but for now set it to a uniform value
            xLims = offset/(2*slope) #identify edges of the scarp
            zNew = np.concatenate((z[x<-xLims], z[x>xLims])) #make copy of all z values outside of the new scarp
            xNew = np.concatenate((x[x<-xLims], x[x>xLims])) #get horizontal positions outside of the scarp
            zNew[xNew<-xLims]+=offset/2 #add half of the offset to the footwall side (just to keep everything centered on zero)
            zNew[xNew>xLims]-=offset/2 #subtract half of the offset to the hangingwall side 
            zTemp = np.interp(x, xNew, zNew) #stitch it together with a linear interpolation over the scarp. This will slightly change the slope depending on dx, but shouldn't create numerical instabilities
            z = zTemp #replace z with zTemp
            eq=0 #set the switch back to zero 
            #pdb.set_trace()
        
        #ESlp = np.random.exponential(1/(E1*slp))
        #E = (E0 + ESlp)*dt*dx
        #E = (E0 + E1*np.abs(slp)**(beta))*dt*dx
        #mat = expMatMake()
        mat = convMatMake()
        Q = np.sum(mat,axis=0)
        dz = -np.gradient(Q)/dx
        #pdb.set_trace()
        #print(np.sum(dz*dx))
        slpTemp = slp+s0
        if np.remainder(t,100)==0:
            print(t)
            xTemp = x/distance
            zTemp = z/distance
            zTemp = (zTemp+s0*xTemp)
            slpTemp = np.gradient(zTemp, xTemp)
            norm = np.trapz(np.abs(slpTemp),xTemp)
            moment1 = np.trapz(xTemp*np.abs(slpTemp)/norm, xTemp)
            moment2 = np.trapz((xTemp-moment1)**2*np.abs(slpTemp)/norm, xTemp)#/np.trapz(slpTemp,x)
            moment3 = np.trapz(((xTemp-moment1)/np.sqrt(moment2))**3*np.abs(slpTemp)/norm, xTemp)#/np.trapz(slpTemp,x)
            Tau = (E1*L0*t)/(distance**2)
            #ax[2,k].plot(t,moment1, 'o', color = clr)
            #ax[3,k].plot(t,moment2, 'o', color = clr)
            ax[2].semilogy(Tau,moment3, 'o', color = clr)
            #if np.remainder(t,100)==0:
             #   ax[1,k].plot(x[np.abs(x)<factor*distance],np.abs(slpTemp[np.abs(x)<factor*distance]), '-', color = clr)
              #  ax[0,k].plot(x[np.abs(x)<factor*distance],z[np.abs(x)<factor*distance], '-', color = clr)

        z+=dz
        t+=dt
        #print(t)
    ax[0].plot(x,z, '-', color = clr)
    ax[1].plot(x,np.abs(slp), '-', color = clr)
    
    
    #ax[1].plot(x,np.abs(slp))

    #beta*=2
    #ax[0,k].plot(x[np.abs(x)<4*distance],z[np.abs(x)<4*distance], '-', color = clr)
    #ax[0,k].set_aspect('equal')
    #ax[0,1].set_aspect('equal')
    #ax[0,2].set_aspect('equal')
    k+=1
    #pdb.set_trace()


#moment1 = np.trapz(x*slpTemp, x)#/np.trapz(slpTemp,x)
#moment2 = np.trapz(x**2*slpTemp, x)#/np.trapz(slpTemp,x)
#moment3 = np.trapz(x**3*slpTemp, x)#/np.trapz(slpTemp,x)
#print(moment1, moment2, moment3)
ax[0].plot(x, fault, '-.', color = 'k')
ax[0].set_ylim(np.min(z)-5, np.max(z)+5)
ax[0].set_ylabel('z')
ax[0].set_xlabel('x')
ax[1].set_ylabel('Slope')
ax[1].set_xlabel('x')
ax[2].set_ylabel('Skewness')
ax[2].set_xlabel('$\\tau = E_1\lambda_0 t/L^2$')
fig.set_tight_layout
fig.set_figwidth(20)

#ax[0].plot(x,z)
#ax[1].plot(x,slpTemp)
plt.show(block = False)
pdb.set_trace()



# ##Fractional Calculus Model ##
# n = 0
# N = 1
# fig, axs = plt.subplots(2,1)
# while n<N:
#     z = np.zeros_like(x)
#     k = np.fft.fftfreq(len(z))

#     count=0
#     num = 100.
#     T = 1000
#     while count<num:
#         offset = np.random.weibull(2.,1)*.5
#         z[np.abs(x)<100.]+=offset
#         ftZ = np.fft.fft(z)
#         ftZT=ftZ*np.exp(K*T*(1j*k)**alpha)
#         z = np.fft.ifft(ftZT)
#         z = np.real(z)
#         T = np.random.weibull(3,1)*100 + 1.
#         #T[T<0] = 1000.
#         print(T)

#         #axs[0].plot(x[int((L-L/2)/dx):int(2*L/dx)], z[int((L-L/2)/dx):int(2*L/dx)])
#         slp = np.gradient(z)/dx

#         #axs[1].plot(x[int((L-L/2)/dx):int(2*L/dx)], slp[int((L-L/2)/dx):int(2*L/dx)])
#         #axs[0].plot(x,z)
#         #axs[1].plot(x,slp)

#         count+=1
#     axs[0].plot(x,z)
#     axs[1].plot(x,slp)
#     n+=1
# plt.show()
