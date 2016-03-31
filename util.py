#!/usr/bin/env python
#############################################################
#EXTERNAL PACKAGES
#############################################################
import numpy as np
from scipy.integrate import odeint
from scipy import constants as const
from math import *
from sys import *

#############################################################
#CONSTANTS
#############################################################

#NUMERICAL CONSTANTS
PI=pi

#PHYSICAL CONSTANTES
MIN=60.0 # s
HOUR=60.0*MIN # s
DAY=24*HOUR # s
YEAR=365.25*DAY # Days
GCONST=const.G # m^3 / kg s^2

#ASTRONOMICAL CONSTANTS
AU=1.496e11 # m
MSUN=1.99e30 # kg

PLANETS=dict(
    Jupiter=dict(
        M= 1.34e27, # kg
        R= 1.47e8, # m
        Prot = 9.4, # s 
        ),
    Saturn=dict(
        M= 1.34e27, # kg
        R= 1.47e8, # m
        Prot = 9.4*HOUR # s 
        ),
    )

#############################################################
#ROUTINES
#############################################################
def Omega(P):
    return 2*PI/P

#############################################################
#TEST MODULE
#############################################################
if __name__=="__main__":
    """
    Normally this routine is used for testing purposes
    """
    print GCONST
    print PLANETS["Jupiter"]["M"]
