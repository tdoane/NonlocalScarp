import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.special import erf
from scipy.optimize import curve_fit
import pdb

def cauchy(x, gamma, alpha):
    #y = alpha*np.exp(-x**2/(2*gamma**2))
    y = alpha/(gamma*(1 + x**2/gamma**2))
    return(y)

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

E0 = 0.000
E1 = 0.005
L0 = .75

k=0
beta = 1.

s0 = 0.1

fig, ax = plt.subplots(1,3)
fig2, ax4 = plt.subplots()
slope = 3.
factor = 30
distance = 0.5

wave = np.fft.fftfreq(len(x))

while k<5:
    distance+= 0.25
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
    tOverD = [0]
    asymm = [0]

    #ax[0, k].plot(x[np.abs(x)<4*distance],z[np.abs(x)<4*distance], '-', color = clrs(0))
    while t<T:
        
        slp = np.gradient(z)/dx
        print(np.min(slp))

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
    
        if np.remainder(t,10)==0:
            print(t)
            xTemp = x#/distance
            zTemp = z#/distance
            zTemp = (zTemp+s0*xTemp)
            slpTemp = np.abs(np.gradient(zTemp, xTemp))
            norm = np.trapz(np.abs(slpTemp),xTemp)
            moment1 = np.trapz(xTemp*np.abs(slpTemp), xTemp)
            moment2 = np.trapz((xTemp-moment1)**2*np.abs(slpTemp), xTemp)#/np.trapz(slpTemp,x)
            moment3 = np.trapz(((xTemp-moment1)/np.sqrt(moment2))**3*np.abs(slpTemp), xTemp)#/np.trapz(slpTemp,x)
            Tau = (E1*L0*t)/(distance**2)

            ftSlp = np.fft.fft(slpTemp)
            reFtSlp = np.real(ftSlp)
            imFtSlp = np.imag(ftSlp)

            ParsIm= np.trapz(np.abs(imFtSlp**2), wave)
            ParsRe = np.trapz(np.abs(reFtSlp**2), wave)
            Pars = np.trapz(np.abs(ftSlp**2), wave)
            ftAsymm = ParsIm/(Pars)
            ax[1].plot(t/(height**2 + distance**2)**(1/2), ftAsymm, 'o', color= clr)

            #ax[2,k].plot(t,moment1, 'o', color = clr)
            #ax[3,k].plot(t,moment2, 'o', color = clr)
            ax[2].plot(t/(height**2 + distance**2)**(1/2), moment3, 'o', color = clr)
            tOverD = np.append(tOverD, t/distance)
            asymm = np.append(asymm, ftAsymm)
            if asymm[1]>asymm[0]:
                asymm  = np.delete(asymm, 0)
                tOverD = np.delete(tOverD, 0)
            #if np.remainder(t,100)==0:
            #   ax[1,k].plot(x[np.abs(x)<factor*distance],np.abs(slpTemp[np.abs(x)<factor*distance]), '-', color = clr)
            #   ax[0,k].plot(x[np.abs(x)<factor*distance],z[np.abs(x)<factor*distance], '-', color = clr)
        z+=dz
        t+=dt
        #print(t)
    ax[0].plot(x,z, '-', color = clr)

    #tOverD = np.delete(tOverD, np.where(asymm<np.max(asymm)))
    #asymm = np.delete(asymm, np.where(asymm<np.max(asymm)))
    #pdb.set_trace()
    #ax[1].plot(x,np.abs(slp), '-', color = clr)
    curve1, curve2 = curve_fit(cauchy, tOverD, asymm, (500, 0.35))
    temp = cauchy(tOverD, curve1[0], curve1[1])
    #ax[1].plot(tOverD, temp, '-k')
    
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
ax[0].set_ylabel('z')
ax[0].set_xlabel('x')
ax[1].set_ylabel('FT Asymmetry')
#ax[1].set_xlabel('$\\tau = E_1\lambda_0 t/L^2$')
ax[1].set_xlabel('$t/h^{1/2}$')
ax[2].set_ylabel('Skewness')
ax[2].set_xlabel('$/h^{1/2}$')
#ax[2].set_xlabel('$\\tau = E_1\lambda_0 t/L^2$')
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
