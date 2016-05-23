#!/usr/bin/env python
# -*- coding:utf-8 -*-
#############################################################
#EXTERNAL PACKAGES
#############################################################
import numpy as np
from scipy.integrate import odeint
from scipy import constants as const
from math import *
from sys import *

#############################################################
#CORE ROUTINES
#############################################################
class dict2obj(object):
    def __init__(self,dic={}):self.__dict__.update(dic)
    def __add__(self,other):
        for attr in other.__dict__.keys():exec("self.%s=other.%s"%(attr,attr))
        return self

#############################################################
#CONSTANTS
#############################################################

#NUMERICAL CONSTANTS
PI=pi
RAD=1/const.degree

#PHYSICAL CONSTANTS
MIN=const.minute # s
HOUR=const.hour # s
DAY=const.day # s
YEAR=const.Julian_year # s
MYEAR=1e6*YEAR # s
GYEAR=1e9*YEAR # s
GCONST=const.G # m^3 / kg s^2

#ASTRONOMICAL CONSTANTS
AU=const.au # m
MSUN=1.99e30 # kg
MENCEL=1.08e20 # kg

PLANETS=dict2obj(dict(
    Jupiter=dict2obj(dict(
                M= 1.898e27, # kg
                R=6.9911e7, # m
                P=29.5*YEAR, # s 
                Prot= 9.4*HOUR, # s 
                alpha=0.126,
                beta=0.020
                )),
    Saturn=dict2obj(dict(
                M=5.683e26, # kg
                R=6.0268e7, # m
                P=10.8*YEAR, # s
                Prot= 10.656*HOUR, # s 
                alpha=0.219,
                beta=0.196
                )),
    Uranus=dict2obj(dict(
                M=86.8e24, # kg
                R=2.5632e7, # m
                P=84*YEAR, # s 
                Prot= 17.24*HOUR, # s 
                alpha=0.30,
                beta=0.093
                )),
    Neptune=dict2obj(dict(
                M=1.024e26, # kg
                R=2.4622e7, # m
                P=164.8*YEAR, # s
                Prot= 16.11*HOUR, # s 
                alpha=0.35,
                beta=0.131
                ))
    ))

#############################################################
#ROUTINES
#############################################################
def omegaAngular(P):
    return 2*PI/P

def semiMajorAxis(P,M):
    a=(GCONST*M*P**2.0/(2.0*PI)**2.0)**(1.0/3.0)
    return a

def meanMotion(a,M):
    n=(GCONST*M/a**3.0)**0.5
    return n

def aRoche(M,densPart=3000,rfac=2.0):
    # Since Roche radius does not depend on R this is a hypotetical one
    R=1.0
    # Planet average density
    densP=M/(4./3*PI*R**3)
    # Roche radius
    ar=rfac*R*(densP/densPart)**(1.0/3.0)
    return ar

def omegaCritic(M,R):
    Oc=np.sqrt(GCONST*M/R**3)
    return Oc

#############################################################
#TEST MODULE
#############################################################
if __name__=="__main__":
    """
    Normally this routine is used for testing purposes
    """
    print GCONST
    print PLANETS.Jupiter.M

    print semiMajorAxis(1*YEAR,MSUN)
    print meanMotion(1*AU,MSUN)*RAD*DAY
    nr=nRoche(PLANETS.Saturn.M,densPart=3000,rfac=2.0)
    print nr

    #print semiMajorAxis(2*PI/nr,1*PLANETS.Jupiter.M)/PLANETS.Jupiter.R

    #print 2*PI/omegaCritic(PLANETS.Jupiter.M,PLANETS.Jupiter.R)/HOUR
    
