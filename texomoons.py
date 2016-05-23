#!/usr/bin/env python
# -*- coding:utf-8 -*-
from util import *

# ############################################################
# CONSTANTS
# ############################################################
# Exponent of radius mass-scaling law
ALPHAR=0.156
# Coefficients of radius time-scaling law
A=3.964
B=-0.064
C=3.364
t0=1e8*YEAR

# Scaling constants for alpha and beta
KP=90.0742985384
BP=4.0
DP=-0.232

#############################################################
# GYRATION RADIUS
#############################################################
GR=0.2 # None

#############################################################
#SPECIFIC TEXOMOONS ROUTINES
#############################################################
def k2Q(alpha,beta,epsilon,**args):
    """
      Source: Mathis, 2015
      alpha = Rc/Rp
      beta = Mc/Mp
      epsilon = Omega/Omega_crit
      args = contains behavior
    """
    #if args["qk2q"]==0:return args["k2q"]

    fac0=alpha**3
    gamma=fac0*(1-beta)/(beta*(1-fac0))

    fac1=alpha**5
    fac2=fac1/(1-fac1)
    fac3=(1-gamma)/gamma*alpha**3

    k2q=100*PI/63*epsilon**2*fac2*(1+fac3)/(1+5./2*fac3)**2
    return k2q

def Mp2Rp(Mp,t,**args):
    """
    args = contains behavior
    """
    if args["qRp"]==0:return args["Rp"]
    Rp=PLANETS.Saturn.R*A*((t+t0)/C)**B

    return Rp

def alpha2beta(Mp,alpha,**args):
    beta=KP*(Mp/PLANETS.Saturn.M)**DP*alpha**BP
    return beta

#############################################################
# DIFFERENTIAL EQUATIONS
#############################################################
def dnmdt(q,t,parameters):
    
    nm=q[0]

    # Primary properties
    Ms=parameters["Ms"]
    Mp=parameters["Mp"]
    alpha=parameters["alpha"]
    Mm=parameters["Mm"]
    ap=parameters["ap"]
    nmp=parameters["nmp"]
    args=parameters["args"]

    # Dynamic parameter
    om=parameters["om"]

    # Secondary properties
    Rp=Mp2Rp(Mp,t,**args)
    epsilon=om/omegaCritic(Mp,Rp)
    beta=alpha2beta(Mp,alpha,**args)
    if args["qk2q"]==0:k2q=args["k2q"]
    else:k2q=k2Q(alpha,beta,epsilon,**args)

    dnmdt=-9./2*k2q*Mm*Rp**5/(GCONST**(5./3)*Mp**(8./3))*\
        nm**(16./3)*np.sign(om-nm)

    return [dnmdt]

def domdt(q,t,parameters):
    
    om=q[0]

    # Primary properties
    Ms=parameters["Ms"]
    Mp=parameters["Mp"]
    alpha=parameters["alpha"]
    Mm=parameters["Mm"]
    ap=parameters["ap"]
    nmp=parameters["nmp"]
    args=parameters["args"]
    
    # Dynamic parameter
    nm=parameters["nm"]

    # Secondary properties
    Rp=Mp2Rp(Mp,t,**args)
    epsilon=om/omegaCritic(Mp,Rp)
    beta=alpha2beta(Mp,alpha,**args)
    if args["qk2q"]==0:k2q=args["k2q"]
    else:k2q=k2Q(alpha,beta,epsilon,**args)

    domdt=-3./2*k2q*Rp**3/(GR*GCONST)*\
        (Mm**2*nm**4*np.sign(om-nm)/Mp**3+\
             (GCONST*Ms)**2*np.sign(om-nmp)/(Mp*ap**6))

    return [domdt]

def dqdt(q,t,parameters):

    nm=q[0]
    om=q[1]

   
    parameters["nm"]=nm
    parameters["om"]=om
    dnmdtp=dnmdt([nm],t,parameters)
    domdtp=domdt([om],t,parameters)

    return dnmdtp+domdtp

#############################################################
#ANALYTICAL SOLUTIONS
#############################################################
def nmt(t,parameters):

    Mm=parameters["Mm"]
    Mp=parameters["Mp"]
    nmini=parameters["nmini"]
    om=parameters["om"]
    alpha=parameters["alpha"]
    args=parameters["args"]

    Rp=Mp2Rp(Mp,t,**args)
    epsilon=om/omegaCritic(Mp,Rp)
    beta=alpha2beta(Mp,alpha)
    if args["qk2q"]==0:k2q=args["k2q"]
    else:k2q=k2Q(alpha,beta,epsilon,**args)
    
    inm=(1/nmini**(13./3)+\
             39./2*k2q*Mm*Rp**5/(GCONST**(5./3)*Mp**(8./3))*t)**(3./13)

    nm=1/inm
    return nm

def omt(t,parameters):

    Ms=parameters["Ms"]
    Mp=parameters["Mp"]
    alpha=parameters["alpha"]
    Mm=parameters["Mm"]
    ap=parameters["ap"]
    args=parameters["args"]

    oini=parameters["oini"]

    Rp=Mp2Rp(Mp,t,**args)
    epsilon=oini/omegaCritic(Mp,Rp)
    beta=alpha2beta(Mp,alpha)
    if args["qk2q"]==0:k2q=args["k2q"]
    else:k2q=k2Q(alpha,beta,epsilon,**args)

    omt=oini-3./2*k2q*GCONST*Rp**3*Ms**2*t/(GR*Mp*ap**6)

    return omt

#############################################################
#TEST
#############################################################
if __name__=="__main__":
    """
    Normally this routine is used for testing purposes
    """
    Mp=PLANETS.Saturn.M
    Rpp=PLANETS.Saturn.R*(Mp/PLANETS.Saturn.M)**ALPHAR
    args=dict(
        qk2q=0,k2q=2.3e-4,
        qRp=0,Rp=Rpp
    )

    print "This is the main module of tExomoons package"
    k2Q=k2Q(0.219,0.126,0.2,**args)
    print "k2Q = ",k2Q

    Rp=Mp2Rp(PLANETS.Saturn.M,0,**args)/PLANETS.Saturn.R
    print "Rp = ",Rp

    """
    #Saturn
    alpha=0.219
    Mp=PLANETS.Saturn.M
    #"""

    #"""
    #Jupiter
    alpha=0.219
    Mp=PLANETS.Saturn.M
    #"""

    """
    #Uranus
    alpha=0.30
    Mp=PLANETS.Uranus.M
    #"""

    """
    #Neptune
    alpha=0.35
    Mp=PLANETS.Neptune.M
    #"""
    
    beta=alpha2beta(Mp,alpha)
    print "Beta = ",beta


