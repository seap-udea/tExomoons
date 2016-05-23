#!/usr/bin/env python
# -*- coding:utf-8 -*-
from texomoons import *
from pylab import*
from math import factorial
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

# ##############################################################
# TIMES AND RADII FOR JUPITER (0.02, 0.1 and 1.0 AU)
# ##############################################################
#Jupiter 1.0 AU
filteredj1AU=np.loadtxt("Points-Jupiter(1AUnc).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tjup1AU=filteredj1AU[:,0]
tjup1AU=tjup1AU.astype(np.float)
Rjup1AU=filteredj1AU[:,1]
Rjup1AU=Rjup1AU.astype(np.float)

filteredj1AUc=np.loadtxt("Points-Jupiter(1AU-25Me).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tjup1AUc=filteredj1AUc[:,0]
tjup1AUc=tjup1AUc.astype(np.float)
Rjup1AUc=filteredj1AUc[:,1]
Rjup1AUc=Rjup1AUc.astype(np.float)

#Jupiter 0.1 AU
filteredj01AU=np.loadtxt("Points-Jupiter(0.1AUnc).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tjup01AU=filteredj01AU[:,0]
tjup01AU=tjup01AU.astype(np.float)
Rjup01AU=filteredj01AU[:,1]
Rjup01AU=Rjup01AU.astype(np.float)

filteredj01AUc=np.loadtxt("Points-Jupiter(0.1AU-25Me).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tjup01AUc=filteredj01AUc[:,0]
tjup01AUc=tjup01AUc.astype(np.float)
Rjup01AUc=filteredj01AUc[:,1]
Rjup01AUc=Rjup01AUc.astype(np.float)

#Jupiter 0.02 AU
filteredj02AU=np.loadtxt("Points-Jupiter(0.02AUnc).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tjup02AU=filteredj02AU[:,0]
tjup02AU=tjup02AU.astype(np.float)
Rjup02AU=filteredj02AU[:,1]
Rjup02AU=Rjup02AU.astype(np.float)

filteredj02AUc=np.loadtxt("Points-Jupiter(0.02AU-25Me).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tjup02AUc=filteredj02AUc[:,0]
tjup02AUc=tjup02AUc.astype(np.float)
Rjup02AUc=filteredj02AUc[:,1]
Rjup02AUc=Rjup02AUc.astype(np.float)

# ##############################################################
# TIMES AND RADII FOR SATURN RESPECTIVELY (0.02, 0.1, 1.0 AU)
# ##############################################################
#Saturn 1.0 AU
filtereds1AU=np.loadtxt("Points-Saturn(1AUnc).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tsat1AU=filtereds1AU[:,0]
tsat1AU=tsat1AU.astype(np.float)
Rsat1AU=filtereds1AU[:,1]
Rsat1AU=Rsat1AU.astype(np.float)*1.16

filtereds1AUc=np.loadtxt("Points-Saturn(1AU-25Me).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tsat1AUc=filtereds1AUc[:,0]
tsat1AUc=tsat1AUc.astype(np.float)
Rsat1AUc=filtereds1AUc[:,1]
Rsat1AUc=Rsat1AUc.astype(np.float)*1.16

#Saturn 0.1 AU
filtereds01AU=np.loadtxt("Points-Saturn(0.1AUnc).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tsat01AU=filtereds01AU[:,0]
tsat01AU=tsat01AU.astype(np.float)
Rsat01AU=filtereds01AU[:,1]
Rsat01AU=Rsat01AU.astype(np.float)*1.16

filtereds01AUc=np.loadtxt("Points-Saturn(0.1AU-25Me).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tsat01AUc=filtereds01AUc[:,0]
tsat01AUc=tsat01AUc.astype(np.float)
Rsat01AUc=filtereds01AUc[:,1]
Rsat01AUc=Rsat01AUc.astype(np.float)*1.16

#Saturn 0.02 AU
filtereds02AU=np.loadtxt("Points-Saturn(0.02AUnc).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tsat02AU=filtereds02AU[:,0]
tsat02AU=tsat02AU.astype(np.float)
Rsat02AU=filtereds02AU[:,1]
Rsat02AU=Rsat02AU.astype(np.float)*1.16

filtereds02AUc=np.loadtxt("Points-Saturn(0.02AU-25Me).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tsat02AUc=filtereds02AUc[:,0]
tsat02AUc=tsat02AUc.astype(np.float)
Rsat02AUc=filtereds02AUc[:,1]
Rsat02AUc=Rsat02AUc.astype(np.float)*1.16

# ##############################################################
# TIMES AND RADII FOR SUPER NEPTUNE (0.02, 0.1, 1.0 AU)
# ##############################################################
#Saturn 1.0 AU
filteredn1AU=np.loadtxt("Points-SuperN(1AUnc).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tnep1AU=filteredn1AU[:,0]
tnep1AU=tnep1AU.astype(np.float)
Rnep1AU=filteredn1AU[:,1]
Rnep1AU=Rnep1AU.astype(np.float)*2.84

filteredn1AUc=np.loadtxt("Points-SuperN(1AU-25Me).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tnep1AUc=filteredn1AUc[:,0]
tnep1AUc=tnep1AUc.astype(np.float)
Rnep1AUc=filteredn1AUc[:,1]
Rnep1AUc=Rnep1AUc.astype(np.float)*2.84

#SNeptune 0.1 AU
filteredn01AU=np.loadtxt("Points-SuperN(0.1AUnc).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tnep01AU=filteredn01AU[:,0]
tnep01AU=tnep01AU.astype(np.float)
Rnep01AU=filteredn01AU[:,1]
Rnep01AU=Rnep01AU.astype(np.float)*2.84

filteredn01AUc=np.loadtxt("Points-SuperN(0.1AU-25Me).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tnep01AUc=filteredn01AUc[:,0]
tnep01AUc=tnep01AUc.astype(np.float)
Rnep01AUc=filteredn01AUc[:,1]
Rnep01AUc=Rnep01AUc.astype(np.float)*2.84

#SNeptune 0.02 AU
filteredn02AU=np.loadtxt("Points-SuperN(0.02AUnc).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tnep02AU=filteredn02AU[:,0]
tnep02AU=tnep02AU.astype(np.float)
Rnep02AU=filteredn02AU[:,1]
Rnep02AU=Rnep02AU.astype(np.float)*2.84

filteredn02AUc=np.loadtxt("Points-SuperN(0.02AU-25Me).csv",\
                        skiprows=1,delimiter=";",dtype='|S9')
tnep02AUc=filteredn02AUc[:,0]
tnep02AUc=tnep02AUc.astype(np.float)
Rnep02AUc=filteredn02AUc[:,1]
Rnep02AUc=Rnep02AUc.astype(np.float)*2.84

# ##############################################################
# CONSTANTS FOR JUPITER, SATURN AND SUPER-NEPTUNE
# ##############################################################
pini=8
pmid=15
pfin=22

#Jupiter constants
alphajup=PLANETS.Jupiter.alpha
Mjup=PLANETS.Jupiter.M
Rjup1AU=Rjup1AU*PLANETS.Jupiter.R
Rjup1AUc=Rjup1AUc*PLANETS.Jupiter.R
Rjup01AU=Rjup01AU*PLANETS.Jupiter.R
Rjup01AUc=Rjup01AUc*PLANETS.Jupiter.R
Rjup02AU=Rjup02AU*PLANETS.Jupiter.R
Rjup02AUc=Rjup02AUc*PLANETS.Jupiter.R
omjup1=2*PI/(pini*HOUR)
omjup2=2*PI/(pmid*HOUR)
omjup3=2*PI/(pfin*HOUR)

#Saturn constants
alphasat=PLANETS.Saturn.alpha
Msat=PLANETS.Saturn.M
Rsat1AU=Rsat1AU*PLANETS.Saturn.R
Rsat1AUc=Rsat1AUc*PLANETS.Saturn.R
Rsat01AU=Rsat01AU*PLANETS.Saturn.R
Rsat01AUc=Rsat01AUc*PLANETS.Saturn.R
Rsat02AU=Rsat02AU*PLANETS.Saturn.R
Rsat02AUc=Rsat02AUc*PLANETS.Saturn.R
omsat1=2*PI/(pini*HOUR)
omsat2=2*PI/(pmid*HOUR)
omsat3=2*PI/(pfin*HOUR)

#SuperN constants
alphanep=PLANETS.Neptune.alpha
Mnep=PLANETS.Neptune.M
Rnep1AU=Rnep1AU*PLANETS.Neptune.R
Rnep1AUc=Rnep1AUc*PLANETS.Neptune.R
Rnep01AU=Rnep01AU*PLANETS.Neptune.R
Rnep01AUc=Rnep01AUc*PLANETS.Neptune.R
Rnep02AU=Rnep02AU*PLANETS.Neptune.R
Rnep02AUc=Rnep02AUc*PLANETS.Neptune.R
omnep1=2*PI/(pini*HOUR)
omnep2=2*PI/(pmid*HOUR)
omnep3=2*PI/(pfin*HOUR)

#print Rsat/PLANETS.Saturn.R

################################################################
# HEAT FUNCTION FOR JUPITER AND SATURN (NO CORE AND WITH CORE)
################################################################
#Jupiter without core
beta=alpha2beta(Mjup,alphajup)
epsilon=omjup1/omegaCritic(Mjup,Rjup1AU)
k2qjup1AU1=k2Q(alphajup,beta,epsilon)
epsilon=omjup2/omegaCritic(Mjup,Rjup1AU)
k2qjup1AU2=k2Q(alphajup,beta,epsilon)
epsilon=omjup3/omegaCritic(Mjup,Rjup1AU)
k2qjup1AU3=k2Q(alphajup,beta,epsilon)

#For Jupiter with 25 Me core
beta=alpha2beta(Mjup,alphajup)
epsilon=omjup1/omegaCritic(Mjup,Rjup1AUc)
k2qjup1AUc1=k2Q(alphajup,beta,epsilon)
epsilon=omjup2/omegaCritic(Mjup,Rjup1AUc)
k2qjup1AUc2=k2Q(alphajup,beta,epsilon)
epsilon=omjup3/omegaCritic(Mjup,Rjup1AUc)
k2qjup1AUc3=k2Q(alphajup,beta,epsilon)

#For Saturn without core
beta=alpha2beta(Msat,alphasat)
epsilon=omsat1/omegaCritic(Msat,Rsat1AU)
k2qsat1AU1=k2Q(alphasat,beta,epsilon)
epsilon=omsat2/omegaCritic(Msat,Rsat1AU)
k2qsat1AU2=k2Q(alphasat,beta,epsilon)
epsilon=omsat3/omegaCritic(Msat,Rsat1AU)
k2qsat1AU3=k2Q(alphasat,beta,epsilon)

#For Saturn with 25 Me core
beta=alpha2beta(Msat,alphasat)
epsilon=omsat1/omegaCritic(Msat,Rsat1AUc)
k2qsat1AUc1=k2Q(alphasat,beta,epsilon)
epsilon=omsat2/omegaCritic(Msat,Rsat1AUc)
k2qsat1AUc2=k2Q(alphasat,beta,epsilon)
epsilon=omsat3/omegaCritic(Msat,Rsat1AUc)
k2qsat1AUc3=k2Q(alphasat,beta,epsilon)

#For SuperN without core
#beta=alpha2beta(Mnep,alphanep)
epsilon=omnep1/omegaCritic(Mnep,Rnep1AU)
k2qnep1AU1=k2Q(alphanep,PLANETS.Neptune.beta,epsilon)
epsilon=omnep2/omegaCritic(Mnep,Rnep1AU)
k2qnep1AU2=k2Q(alphanep,PLANETS.Neptune.beta,epsilon)
epsilon=omnep3/omegaCritic(Mnep,Rnep1AU)
k2qnep1AU3=k2Q(alphanep,PLANETS.Neptune.beta,epsilon)

#For SuperN with 25 Me core
#beta=alpha2beta(Mnep,alphanep)
epsilon=omnep1/omegaCritic(Mnep,Rnep1AUc)
k2qnep1AUc1=k2Q(alphanep,PLANETS.Neptune.beta,epsilon)
epsilon=omnep2/omegaCritic(Mnep,Rnep1AUc)
k2qnep1AUc2=k2Q(alphanep,PLANETS.Neptune.beta,epsilon)
epsilon=omnep3/omegaCritic(Mnep,Rnep1AUc)
k2qnep1AUc3=k2Q(alphanep,PLANETS.Neptune.beta,epsilon)

##############################################################
#FIGURES OF PLANETARY HEAT FUNCTION 
##############################################################

#Jupiter without core
smoothj1AU1=savitzky_golay(k2qjup1AU1, 51, 4) # window size 51, polynomial order 3
smoothj1AU2=savitzky_golay(k2qjup1AU2, 51, 4) # window size 51, polynomial order 3
smoothj1AU3=savitzky_golay(k2qjup1AU3, 51, 4) # window size 51, polynomial order 3

#Jupiter with core
smoothj1AUc1=savitzky_golay(k2qjup1AUc1, 51, 4) # window size 51, polynomial order 3
smoothj1AUc2=savitzky_golay(k2qjup1AUc2, 51, 4) # window size 51, polynomial order 3
smoothj1AUc3=savitzky_golay(k2qjup1AUc3, 51, 4) # window size 51, polynomial order 3

#Saturn without core
smooths1AU1=savitzky_golay(k2qsat1AU1, 51, 4) # window size 51, polynomial order 3
smooths1AU2=savitzky_golay(k2qsat1AU2, 51, 4) # window size 51, polynomial order 3
smooths1AU3=savitzky_golay(k2qsat1AU3, 51, 4) # window size 51, polynomial order 3

#For Saturn with core
smooths1AUc1=savitzky_golay(k2qsat1AUc1, 51, 4) # window size 51, polynomial order 3
smooths1AUc2=savitzky_golay(k2qsat1AUc2, 51, 4) # window size 51, polynomial order 3
smooths1AUc3=savitzky_golay(k2qsat1AUc3, 51, 4) # window size 51, polynomial order 3

#For Super-Neptune without core
smoothn1AU1=savitzky_golay(k2qnep1AU1, 51, 4) # window size 51, polynomial order 3
smoothn1AU2=savitzky_golay(k2qnep1AU2, 51, 4) # window size 51, polynomial order 3
smoothn1AU3=savitzky_golay(k2qnep1AU3, 51, 4) # window size 51, polynomial order 3

#For Super-Neptune with core
smoothn1AUc1=savitzky_golay(k2qnep1AUc1, 51, 4) # window size 51, polynomial order 3
smoothn1AUc2=savitzky_golay(k2qnep1AUc2, 51, 4) # window size 51, polynomial order 3
smoothn1AUc3=savitzky_golay(k2qnep1AUc3, 51, 4) # window size 51, polynomial order 3

fig=plt.figure()
gs = gridspec.GridSpec(1, 2)
gs.update(wspace=0.8,top=0.9)
#fig.suptitle("$-$ No core, -- 25 $M_{\oplus}$ core",fontsize=14,y=1.001)
#ax1=plt.subplot(121)
ax1=plt.subplot(gs[0, :1], )
#ax=subplot2grid((2,2),(0,0))
ax1.plot(tjup1AU*YEAR/GYEAR,smoothj1AU1,'r-')
ax1.plot(tjup1AU*YEAR/GYEAR,smoothj1AU2,'g-')
ax1.plot(tjup1AU*YEAR/GYEAR,smoothj1AU3,'b-')
ax1.plot(tjup1AUc*YEAR/GYEAR,smoothj1AUc1,'r--')
ax1.plot(tjup1AUc*YEAR/GYEAR,smoothj1AUc2,'g--')
ax1.plot(tjup1AUc*YEAR/GYEAR,smoothj1AUc3,'b--')
#ax1.legend(loc='best',fontsize=12,frameon=True)
ax1.set_title("Jupiter (0.1 AU)",fontsize=12)
ax1.set_ylabel(r'$k_{2}/Q_{p}$',fontsize=15)
ax1.set_xlabel('$\mathrm{Time\;[Gyrs]}$',fontsize=13)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlim((0,4.6))
plt.yticks(fontsize=11)
plt.text(1,0.73e-4, '$-$ No core\n -- 25 $M_{\oplus}$ core',\
             ha='center',va='center',size=12, alpha=.8,bbox=dict(facecolor='none',edgecolor='black',pad=4.0))

#ax2=plt.subplot(122)
#ax2=subplot2grid((2,2),(0,1))
ax2=plt.subplot(gs[0, 1:])
ax2.plot(tsat1AU*YEAR/GYEAR,smooths1AU1,'r-',label='$P_{\mathrm{rot}}=8 h$')
ax2.plot(tsat1AU*YEAR/GYEAR,smooths1AU2,'g-',label='$P_{\mathrm{rot}}=15 h$')
ax2.plot(tsat1AU*YEAR/GYEAR,smooths1AU3,'b-',label='$P_{\mathrm{rot}}=22 h$')
ax2.plot(tsat1AUc*YEAR/GYEAR,smooths1AUc1,'r--')
ax2.plot(tsat1AUc*YEAR/GYEAR,smooths1AUc2,'g--')
ax2.plot(tsat1AUc*YEAR/GYEAR,smooths1AUc3,'b--')
ax2.legend(loc='best',fontsize=12,frameon=True)
ax2.set_title("Saturn (0.1 AU)",fontsize=12)
ax2.set_ylabel(r'$k_{2}/Q_{p}$',fontsize=15)
ax2.set_xlabel('$\mathrm{Time\;[Gyrs]}$',fontsize=13)
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xlim((0,4.6))
plt.yticks(fontsize=11)

'''#ax=plt.subplot(213)
#ax=subplot2grid((2,2),(1,0))
ax3=plt.subplot(gs[1, 1:3])
ax3.plot(tnep1AU*YEAR/GYEAR,smoothn1AU1,'k-')
ax3.plot(tnep1AU*YEAR/GYEAR,smoothn1AU2,'r-')
ax3.plot(tnep1AU*YEAR/GYEAR,smoothn1AU3,'b-')
ax3.plot(tnep1AUc*YEAR/GYEAR,smoothn1AUc1,'k--')
ax3.plot(tnep1AUc*YEAR/GYEAR,smoothn1AUc2,'r--')
ax3.plot(tnep1AUc*YEAR/GYEAR,smoothn1AUc3,'b--')
ax3.legend(loc='best',fontsize=12,frameon=True)
ax3.set_title("Super-Neptune")
ax3.set_ylabel(r'$k_{2}/Q_{p}$',fontsize=15)
ax3.set_xlabel('$\mathrm{Time\;[Gyrs]}$',fontsize=15)
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.set_xlim((0,4.6))
plt.yticks(fontsize=11)
'''
gs.tight_layout(fig, rect=[None, 0, 1, 1], h_pad=0.1)
plt.savefig("variationk2Q-1AU.png")

##############################################################
#FIGURES OF PLANETARY RADIUS EVOLUTION
##############################################################

#JUPITER

radsmoothj1AU=savitzky_golay(Rjup1AU, 51, 4) # window size 51, polynomial order 3
radsmoothj1AUc=savitzky_golay(Rjup1AUc, 51, 4) # window size 51, polynomial order 3

radsmoothj01AU=savitzky_golay(Rjup01AU, 51, 4) # window size 51, polynomial order 3
radsmoothj01AUc=savitzky_golay(Rjup01AUc,51,4) # window size 51, polynomial order 3

radsmoothj02AU=savitzky_golay(Rjup02AU, 51, 4) # window size 51, polynomial order 3
radsmoothj02AUc=savitzky_golay(Rjup02AUc, 51, 4) # window size 51, polynomial order 3

fig=plt.figure()
gs = gridspec.GridSpec(2, 4)
gs.update(wspace=0.8,top=0.8)
#fig.suptitle("$-$ No core, -- 25 $M_{\oplus}$ core",fontsize=14,y=1.001)
#ax=plt.subplot(221)
ax1=plt.subplot(gs[0, :2], )
#ax=subplot2grid((2,2),(0,0))
ax1.plot(tjup1AU*YEAR/GYEAR,radsmoothj1AU/PLANETS.Jupiter.R,'r-')
ax1.plot(tjup1AUc*YEAR/GYEAR,radsmoothj1AUc/PLANETS.Jupiter.R,'r--')
ax1.plot(tjup01AU*YEAR/GYEAR,radsmoothj01AU/PLANETS.Jupiter.R,'g-')
ax1.plot(tjup01AUc*YEAR/GYEAR,radsmoothj01AUc/PLANETS.Jupiter.R,'g--')
ax1.plot(tjup02AU*YEAR/GYEAR,radsmoothj02AU/PLANETS.Jupiter.R,'b-')
ax1.plot(tjup02AUc*YEAR/GYEAR,radsmoothj02AUc/PLANETS.Jupiter.R,'b--')
#ax1.legend(loc='best',fontsize=12,frameon=True)
ax1.set_title("Jupiter",fontsize=12)
ax1.set_ylabel(r'$R_{\mathrm{p}}\;(R_{\mathrm{j}})$',fontsize=13)
ax1.set_xlabel('$\mathrm{Time\;(Gyrs)}$',fontsize=13)
ax1.set_xlim((0,6))
plt.yticks(fontsize=11)
plt.text(0.5,1.738, 'I', ha='center', va='center',size=12,color='r',alpha=.8)
plt.text(4.45,1.658, '$-$ No core\n -- 25 $M_{\oplus}$ core',\
             ha='center',va='center',size=12, alpha=.8,bbox=dict(facecolor='none',edgecolor='black',pad=4.0))

#SATURN

radsmooths1AU=savitzky_golay(Rsat1AU, 51, 4) # window size 51, polynomial order 3
radsmooths1AUc=savitzky_golay(Rsat1AUc, 51, 4) # window size 51, polynomial order 3

radsmooths01AU=savitzky_golay(Rsat01AU, 51, 4) # window size 51, polynomial order 3
radsmooths01AUc=savitzky_golay(Rsat01AUc,51,4) # window size 51, polynomial order 3

radsmooths02AU=savitzky_golay(Rsat02AU, 51, 4) # window size 51, polynomial order 3
radsmooths02AUc=savitzky_golay(Rsat02AUc, 51, 4) # window size 51, polynomial order 3

#ax=plt.subplot(222)
#ax=subplot2grid((2,2),(0,1))
ax2=plt.subplot(gs[0, 2:])
ax2.plot(tsat1AU*YEAR/GYEAR,radsmooths1AU/PLANETS.Saturn.R,'r-',\
             label=r'1.0 AU')
ax2.plot(tsat1AUc*YEAR/GYEAR,radsmooths1AUc/PLANETS.Saturn.R,'r--')
ax2.plot(tsat01AU*YEAR/GYEAR,radsmooths01AU/PLANETS.Saturn.R,'g-',\
             label=r'0.1 AU')
ax2.plot(tsat01AUc*YEAR/GYEAR,radsmooths01AUc/PLANETS.Saturn.R,'g--')
ax2.plot(tsat02AU*YEAR/GYEAR,radsmooths02AU/PLANETS.Saturn.R,'b-',\
             label=r'0.02 AU')
ax2.plot(tsat02AUc*YEAR/GYEAR,radsmooths02AUc/PLANETS.Saturn.R,'b--')
ax2.legend(loc='best',fontsize=12,frameon=True)
ax2.set_title("Saturn",fontsize=12)
ax2.set_ylabel(r'$R_{\mathrm{p}}\;(R_{\mathrm{Sat}})$',fontsize=13)
ax2.set_xlabel('$\mathrm{Time\;(Gyrs)}$',fontsize=13)
ax2.set_xlim((0,6))
plt.yticks(fontsize=11)
plt.text(0.5,2.838, 'II', ha='center', va='center',size=12,color='r',alpha=.8)
#plt.text(0.7,2.9, 'Saturn', ha='center', va='center',
        #size=12, alpha=.8)

#NEPTUNE

radsmoothn1AU=savitzky_golay(Rnep1AU, 51, 4) # window size 51, polynomial order 3
radsmoothn1AUc=savitzky_golay(Rnep1AUc, 51, 4) # window size 51, polynomial order 3

radsmoothn01AU=savitzky_golay(Rnep01AU, 51, 4) # window size 51, polynomial order 3
radsmoothn01AUc=savitzky_golay(Rnep01AUc,51,4) # window size 51, polynomial order 3

radsmoothn02AU=savitzky_golay(Rnep02AU, 51, 4) # window size 51, polynomial order 3
radsmoothn02AUc=savitzky_golay(Rnep02AUc, 51, 4) # window size 51, polynomial order 3

#ax=plt.subplot(213)
#ax=subplot2grid((2,2),(1,0))
ax3=plt.subplot(gs[1, 1:3])
ax3.plot(tnep1AU*YEAR/GYEAR,radsmoothn1AU/PLANETS.Neptune.R,'r-')
ax3.plot(tnep1AUc*YEAR/GYEAR,radsmoothn1AUc/PLANETS.Neptune.R,'r--')
ax3.plot(tnep01AU*YEAR/GYEAR,radsmoothn01AU/PLANETS.Neptune.R,'g-')
ax3.plot(tnep01AUc*YEAR/GYEAR,radsmoothn01AUc/PLANETS.Neptune.R,'g--')
ax3.plot(tnep02AU*YEAR/GYEAR,radsmoothn02AU/PLANETS.Neptune.R,'b-')
ax3.plot(tnep02AUc*YEAR/GYEAR,radsmoothn02AUc/PLANETS.Neptune.R,'b--')
#ax3.legend(loc='best',fontsize=12,frameon=True)
ax3.set_title("Super Neptune ($\sim2\,M_{\mathrm{Nep}}$)",fontsize=12)
ax3.set_ylabel(r'$R_{\mathrm{p}}\;(R_{\mathrm{Nep}})$',fontsize=13)
ax3.set_xlabel('$\mathrm{Time\;(Gyrs)}$',fontsize=13)
ax3.set_xlim((0,6))
plt.yticks(fontsize=11)
plt.text(0.5,6.55, 'III', ha='center', va='center',size=12,color='r',alpha=.8)
gs.tight_layout(fig, rect=[None, 0, 1, 1], h_pad=0.1)
plt.savefig("Radius-evolution-threeplanets.png")


'''
plt.figure()
plt.plot(tjup02AUc*YEAR/GYEAR,radsmoothj02AUc/PLANETS.Jupiter.R,'b--')
plt.show()
'''
