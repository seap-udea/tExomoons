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

#PHYSICAL CONSTANTS
MIN=60.0 # s
HOUR=60.0*MIN # s
DAY=24*HOUR # s
YEAR=365.25*DAY # s
MYEAR=YEAR**8.0 # s
GYEAR=YEAR**9.0 # s
GCONST=const.G # m^3 / kg s^2
GR=0.2 # None
DensPart=3000 # kg / m^3

#ASTRONOMICAL CONSTANTS
AU=1.496e11 # m
MSUN=1.99e30 # kg
MENCEL=1.08e20 # kg
MSAT=5.683e26 # kg
MJUP=1.898e27 # kg
MNEP=1.024e26 # kg
MURAN=86.8e24 # kg
RSAT=60268000 # m
RJUP=69911000 # m
RNEP=24622000 # m
RURAN=25632000 # m
PSAT=10759.22*DAY # s
PJUP=4332.6*DAY # s
PNEP=60189*DAY # s
PURAN=30685.4*DAY # s

#INTEGRATION CONSTANTS
Nts=100*GYEAR
dt=0.001*GYEAR

#PARAMETERS
PR=15.0*HOUR # s
PO=70.0*DAY # s
MS=1.0*MSUN # kg
MP=1.0*MSAT # kg
MM=10.0*MENCEL # kg

PLANETS=dict(
    Jupiter=dict(
        M= MJUP, # kg
        R= RJUP, # m
        Porb= PJUP, # s
        Prot= 9.4*HOUR, # s 
        alpha=0.126,
        beta=0.020
        ),
    Saturn=dict(
        M= MSAT, # kg
        R= RSAT, # m
        Porb= PSAT, # s
        Prot= 10.656*HOUR, # s 
        alpha=0.219,
        beta=0.196
        ),
    Uranus=dict(
        M= MURAN, # kg
        R= RURAN, # m
        Porb=PURAN, # s
        Prot= 17.24*HOUR, # s 
        alpha=0.30,
        beta=0.093
        ),
    Neptune=dict(
        M= MNEP, # kg
        R= RNEP, # m
        Porb= PNEP, # s
        Prot= 16.11*HOUR, # s 
        alpha=0.35,
        beta=0.131
        )
    )

#############################################################
#ROUTINES
#############################################################
def Omega(PR):
    return 2*PI/PR

def Semimaj(PO):
    return ((GCONS*MS*PO**2.0)/(2.0*PI)**2.0)**(1.0/3.0)

def nroche(MP):
    DensP=(3.0*MP)/(4*PI)
    B1=(DensP/DensPart)**(1.0/3.0)
    ar=2.0*B1
    return ((GCONS*MP)/ar**3.0)**0.5
    

#############################################################
#TEST MODULE
#############################################################
if __name__=="__main__":
    """
    Normally this routine is used for testing purposes
    """
    print GCONST
    print PLANETS["Jupiter"]["M"]
