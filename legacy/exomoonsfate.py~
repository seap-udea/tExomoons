import numpy as np
from matplotlib.pyplot import*
from scipy.integrate import odeint
from math import*
from sys import*

###################CONSTANTS#####################
N=3.0

#INPUT PARAMETERS
amu=3.5;"Initial moon position (Rp)";
tsu=80.0;"Planet orbital period (Days)";
alph=0.219;"alpha of Saturn";
#bet=0.196
#bet=0.263
A=3.964
B=-0.064
C=3.364

#######COEFICIENTS FROM THE SCALE-LAW FOR THE PLANET'S TIDAL QUALITY FACTOR
kpp=259.179864469
bp=4.0
delta=-0.231999999999998

#############################
#FOR THE MASSES CONTOURS
#############################
Minip=0.5;"Minimum planet mass (Msat)";
Maxip=3.0;"Maximum planet mass (Msat)";

Minim=1.0;"Minimum moon mass (Mencel)";
Maxim=12.0;"Maximum moon mass (Mencel)";

Mps=np.linspace(Minip,Maxip,N);"Planet Mass (Msat)";
Mms=np.linspace(Minim,Maxim,N);"Moon Mass (Mencel)";
MMS,MPS=np.meshgrid(Mms,Mps)


############################
#INITIAL CONDITIONS
############################
tini=0.0
#################################################
#SWITCH THE EFFECT OF PLANETARY RADIUS EVOLUTION
################################################

#band=0;"The radius of the planet evolves while the moon changes its orbit";
band=1;"The radius of the planet is constant throghout the round-trip";

########################################################
#TIME VECTOR DEPENDING ON EVOLUTION OF PLANETARY RADIUS
########################################################

ts=np.arange(tini,Nts*dt,dt)

############################
#RUNGE-KUTTA 4
############################

def eom(q,t,parameters):
    band=parameters['band']
    Oini=parameters['Oini']
    Gyear=parameters['Gyear']
    Msat=parameters['Msat']
    G=parameters['G']
    GR=parameters['GR']
    t0=parameters['t0']
    A=parameters['A']
    B=parameters['B']
    C=parameters['C']
    Mp=parameters['Mp']
    Mm=parameters['Mm']
    am=parameters['am']
    ap=parameters['ap']
    alph=parameters['alph']
    delta=parameters['delta']
    omega=q[0]*0+Oini*1
    motion=q[1]
    dtheta=omega
    dmotion=motion

    # Depending on Mp and Mm there are several quantities
    if band==0:
        Rp=A*((t+t0)/C)**B
        
    elif band==1:
        Rp=(Mp/Msat)**0.156

    bet=kpp*alph**bp*(Mp*95.159)**delta;"Beta of the planet";
    Y=(alph**3.0*(1.0-bet))/(bet*(1.0-alph**3.0))
    T1=100.0*pi/63.0
    T2=alph**5.0/(1.0-alph**5.0)
    T3=1.0+((1.0-Y)/Y)*alph**3.0
    T4=(1.0+(2.5*(1.0-Y)/Y)*alph**3.0)**(-2.0)
    Cons=T1*T2*T3*T4
    
    eps2=omega**2.0/(G*Mp/Rp**3.0)
    KQ=eps2*Cons
    
    Num=3.0*KQ*G*Rp**3.0*Ms**2.0
    Den=2.0*GR*Mp*ap**6.0
    TORS=Num/Den
    
    NUME=9*KQ*Mm*Rp**5.0
    DENO=2*G**(5.0/3.0)*Mp**(8.0/3.0)
    TORP=NUME/DENO
    
    Numer=3.0*KQ*Rp**3.0*Mm**2.0
    Denom=2.0*GR*G*Mp**3.0
    TORM=Numer/Denom

    DOmegap=-TORS-TORM*motion**4.0
    Dmotion=-TORP*motion**(16.0/3.0)

    return np.array([DOmegap,Dmotion])

doms=np.zeros((N,N))
dnms=np.zeros((N,N))
tsync=np.zeros((N,N))
for i in np.arange(N):
    for j in np.arange(N):
        #For each i,j I have a Mp and an Mm
        Mpij=Mp[i,j]
        Mmij=Mm[i,j]

        #There is a variety of quantities depending on Mp and Ms

        Rpp=(Mpij/Msat)**0.156
        am=3.5*Rpp
        nini=((G*Mpij)/am**3.0)**(0.5);"initial orbital lunar mean motion";

        # This is the solution given an nini, Mp, Mm, and a set of constants

        e=[Oini,nini]
        pars=dict(Oini=Oini,Mp=Mpij,Mm=Mmij,ap=ap,am=am,A=A,B=B,C=C,GR=GR,G=G,Gyear=Gyear,Msat=Msat,t0=t0,alph=alph,delta=delta,band=band)
        solution=odeint(eom,e,ts*Gyear,args=(pars,))
        doms=solution[:,0]
        dnms=solution[:,1]
            
        print doms[0:10]
        print dnms[0:10]

exit(0)
############################################################
#WE GET THE SYNCHRONIZATION TIME FOR THE PLANET-MOON SYSTEM
############################################################

tsync=[]
for j in np.arange(len(doms)):
    if (doms[j]<=dnms[j]):
        tsync+=[ts[j]]
tsync=np.array(tsync)

#Tsynch=tsync.min();"synchronous time for band=0";

print tsync
exit(0)
###############################
#PLOT
###############################

figure(figsize=(8.0,8.0))

plot(ts,doms,label=r'$\Omega_{\mathrm{p}}(t)$')
plot(ts,dnms,label='$n_{\mathrm{m}}^{-}(t)$')
#title(r"$k_{\mathrm{p}}/Q(R_{\mathrm{p}}(t)\,,\,\Omega_{\mathrm{p,cte}})$ , $P_{\mathrm{orb}}=%.1f$ d , $\alpha=%.3f$ , $\beta=%.3f$"%(tsu,alph,bet))
#title(r"$k_{\mathrm{p}}/Q(R_{\mathrm{p,cte}}\,,\,\Omega_{\mathrm{p,cte}})$ , $P_{\mathrm{orb}}=%.1f$ d , $\alpha=%.3f$ , $\beta=%.3f$"%(tsu,alph,bet))
legend(loc='best',frameon=False)
xlabel("Time [Gyr]")
ylabel(r"$\Omega_{\mathrm{p}}$ , $n_{\mathrm{m}}$ [yr$^{-1}$]")
#savefig("evolutionrate-radvaromegaconst.png")
savefig("evolutionrate-radconstomegaconst.png")


########################################################
#WE GET THE DECAY TIME FOR THE PLANET-MOON SYSTEM
########################################################
#Nms=1028;"Rp=cte,Om=cte";
#Nms=154.1;"Rp=var,Om=cte";
Nms=90;"Rp=var,Om=cte";
#Nms=4000;"Rp=var,Om=int";
#Nms=5000;"Rp=cte,Om=int";

dtt=0.001;"time-step";

############################
#INITIAL CONDITIONS
############################

Oinitial=doms[Tsynch/dt]
ninitial=dnms[Tsynch/dt]

############################
#RUNGE-KUTTA 4
############################

def dec(q,t,parameters):
    Oi=parameters['Oi']
    Mp=parameters['Mp']
    band=parameters['band']
    omega=q[0]*0+Oi*1
    motion=q[1]
    dtheta=omega
    dmotion=motion
    DOmegap=-K2Q(band,omega,t)[1]+K2Q(band,omega,t)[0]*motion**4.0
    Dmotion=+K2Q(band,omega,t)[2]*motion**(16.0/3.0)
    return np.array([DOmegap,Dmotion])

############################
#SOLUTION
############################

time=np.arange(tini,Nms*dtt,dtt)
time=np.array(time)

IC=[Oinitial,ninitial]
pars=dict(Oi=Oi,Mp=Mp,band=band)
solve,info=odeint(dec,IC,time*Gyear,args=(pars,),full_output=True)

oms=solve[:,0]
nms=solve[:,1]

print npl

##############################################################
#WE CALCULATE THE TIME EVOLUTION OF THE LUNAR POSITION (DECAY)
##############################################################
adecay=((G*Mp)/nms**2.0)**(1.0/3.0)


figure(figsize=(7.5,7.5))

plot(time,nms,'-g',label='Lunar orbital decay until $a_{\mathrm{Roche}}$')
#title(r"$k_{\mathrm{p}}/Q(R_{\mathrm{p,cte}}\,,\,\Omega_{\mathrm{p}}(t))$ , $P_{\mathrm{orb}}=%.1f$ d , $\alpha=%.3f$ , $\beta=%.3f$"%(tsu,alph,bet))
#title(r"$k_{\mathrm{p}}/Q(R_{\mathrm{p}}(t)\,,\,\Omega_{\mathrm{p}}(t))$ , $P_{\mathrm{orb}}=%.1f$ d , $\alpha=%.3f$ , $\beta=%.3f$"%(tsu,alph,bet))
title(r"$k_{\mathrm{p}}/Q(R_{\mathrm{p}}(t)\,,\,\Omega_{\mathrm{p,cte}})$ , $P_{\mathrm{orb}}=%.1f$ d , $\alpha=%.3f$ , $\beta=%.3f$"%(tsu,alph,bet))
#title(r"$k_{\mathrm{p}}/Q(R_{\mathrm{p,cte}}\,,\,\Omega_{\mathrm{p,cte}})$ , $P_{\mathrm{orb}}=%.1f$ d , $\alpha=%.3f$ , $\beta=%.3f$"%(tsu,alph,bet))
axhline(nroche,ls='--',color='r',linewidth=1.5,label='$n_{\mathrm{Roche}}$')
legend(loc='best',frameon=False)
xlabel("Time [Gyr]")
ylabel(r"$n_{\mathrm{m}}$ [yr$^{-1}$]")
ylim((ninitial,11000))
#savefig("motionevolution-radconstomegaint.png")
#savefig("motionevolution-radvaromegaint.png")
savefig("motionevolution-radvaromegaconst.png")
#savefig("motionevolution-radconstomegaconst.png")

figure(figsize=(7.5,7.5))

plot(vec,async,'k-',label='Lunar orbital decay until $a_{\mathrm{Roche}}$')
plot(time+Tsynch,adecay,'k-')
#title(r"$k_{\mathrm{p}}/Q(R_{\mathrm{p,cte}}\,,\,\Omega_{\mathrm{p}}(t))$ , $P_{\mathrm{orb}}=%.1f$ d , $\alpha=%.3f$ , $\beta=%.3f$"%(tsu,alph,bet))
#title(r"$k_{\mathrm{p}}/Q(R_{\mathrm{p}}(t)\,,\,\Omega_{\mathrm{p}}(t))$ , $P_{\mathrm{orb}}=%.1f$ d , $\alpha=%.3f$ , $\beta=%.3f$"%(tsu,alph,bet))
title(r"$k_{\mathrm{p}}/Q(R_{\mathrm{p}}(t)\,,\,\Omega_{\mathrm{p,cte}})$ , $P_{\mathrm{orb}}=%.1f$ d , $\alpha=%.3f$ , $\beta=%.3f$"%(tsu,alph,bet))
#title(r"$k_{\mathrm{p}}/Q(R_{\mathrm{p,cte}}\,,\,\Omega_{\mathrm{p,cte}})$ , $P_{\mathrm{orb}}=%.1f$ d , $\alpha=%.3f$ , $\beta=%.3f$"%(tsu,alph,bet))
axhline(ar,ls='--',color='r',linewidth=1.5,label='$a_{\mathrm{Roche}}$')
legend(loc='best',frameon=False)
xlabel("Time [Gyr]",fontsize=17)
ylabel("$a_{\mathrm{m}}/R_{\mathrm{p}}$",fontsize=17)
#savefig("orbitaldecay-radconstomegaint.png")
#savefig("orbitaldecay-radvaromegaint.png")
savefig("orbitaldecay-radvaromegaconst.png")
#savefig("orbitaldecay-radconstomegaconst.png")

tdec=[]
for i in np.arange(len(nms)):
    if nms[i]>=nroche:
        tdec+=[time[i]]
tdec=np.array(tdec)

Tdecay=tdec.min();"decay time0";

Trt=Tsynch+Tdecay

print  "Tdecay= %.3f Gyr, Round-trip= %.3f Gyr "%(Tdecay,Trt)


#THESE ARE ROUTINES FOR THE CONTOURS WITH ANALYTICAL FORMULAE
#WE SET UP THE NUMBER OF ITERATIONS AND THE INITIAL GUESS. WE CALL
#THE FUNCTION "tsync" WITH "tini" AND "numiter" (WE MUST DO THIS FOR RUN
#THE PROGRAM)

tini=1e5
numiter=10

#WE DEFINE HERE THE FORM OF THE ARRANGED EQUATION SOLVING SUCH THAT 't=g(t,Mp,Mm,ap,am)'

def g(t,P,Ms,Mp,Mm,am,A):
    np=(2*pi)/P
    ap=((G*Ms)/np**(2.0))**(1.0/3.0)
    Rp=(Mp/Msat)**(0.156)
    eps2=(Oi/(G*Mp/Rp**3.0)**(0.5))**2.0
    B=kpp*A**bp*(Mp*95.159)**delta;"Beta of the planet";
    #B=0.196
    Y=(A**3.0*(1-B))/(B*(1-A**3.0));"Gamma of the planet";
    T1=100*pi/63
    T2=A**5.0/(1-A**5.0)
    T3=1+((1-Y)/Y)*A**3.0
    T4=(1+(2.5*(1-Y)/Y)*A**3.0)**(-2.0)
    KQ=T1*eps2*T2*T3*T4
    #print (KQ)
    nmi=((G*Mp)/am**3.0)**(0.5)
    A1=(nmi**(-13.0/3.0)+(19.5*KQ*Mm*(Rp)**(5.0)*t)/(G**(5.0/3.0)*(Mp)**(8.0/3.0)))
    A2=Oi-(A1)**(-3.0/13.0)
    A3=A2*((2*GR*(Mp))/(3*KQ*G*Rp**(3.0)))
    A4=(Ms**(2.0)/ap**(6.0)+(Mm)**(2.0)/am**(6.0))
    Gfun=A3*A4**(-1.0)
    return Gfun

#THIS IS THE FIXED-POINT FUNCTION FOR DOING THE "k" ITERATIONS TO
#REFINE THE SEED'S VALUE (tini).

def syncTime(t,k,P,Ms,Mp,Mm,am,A):
    for i in np.arange(k):
        t=g(t,P,Ms,Mp,Mm,am,A)
    return t

#THIS IS GOING TO BE THE DATA FOR DOING THE CALCULUS OF nsync AND async

def decayTime(P,Ms,Mp,Mm,am,tsync,A):
    np=(2*pi)/P
    ap=((G*Ms)/np**(2.0))**(1.0/3.0)
    Rp=(Mp/Msat)**(0.156)
    eps2=(Oi**2.0/(G*Mp/Rp**3.0))
    #B=0.196
    B=kpp*A**bp*(Mp*95.159)**delta;"Beta of the planet";
    Y=(A**3.0*(1-B))/(B*(1-A**3.0));"Gamma of the planet";
    T1=100*pi/63
    T2=A**5.0/(1-A**5.0)
    T3=1+((1-Y)/Y)*A**3.0
    T4=(1+(2.5*(1-Y)/Y)*A**3.0)**(-2.0)
    KQ=T1*eps2*T2*T3*T4
    nmi=((G*(Mp))/am**(3.0))**(0.5)
    NMI=nmi**(-13.0/3.0)
    
    #CALCULUS OF THE ROCHE RADIUS
    DensP=(3.0*(Mp))/(4*pi*Rp**(3.0))
    B1=(DensP/DensPart)**(1.0/3.0)
    ar=2.0*B1*Rp
    nroche=((G*Mp)/ar**(3.0))**0.5
    #print ("Roche Radius= %.3f Rsat"%(ar))
    
    NUM=(19.5*KQ*(Mm)*Rp**(5.0)*tsync)
    DEN=(G**(5.0/3.0)*(Mp)**(8.0/3.0))
    nsync=(NMI+NUM/DEN)**(-3.0/13.0)
    async=((G*(Mp))/nsync**(2.0))**(1.0/3.0)

    #NEXT, WE CALCULATE THE DECAY TIME FROM THE SYNCHRONOUS RADIUS TO THE
    #ROCHE RADIUS
    TER1=(19.5*KQ*(Mm)*(Rp)**(5.0))
    TER2=(G**(5.0/3.0)*(Mp)**(8.0/3.0))
    TER3=TER2/TER1
    TER4=(nsync**(-13.0/3.0)-nroche**(-13.0/3.0))
    Td=TER3*TER4
    return Td



###############CONTOURS FOR THE PLANETARY MASSES##################
##################################################################

#FIXED PARAMETERS
Msu=1.0;"Stellar mass (Msun)";
amu=4.0;"Initial moon position (Rp)";
tsu=70.0;"Planet orbital period (Days)";
A=0.18;"Alpha of the planet";

#FOR THE MASSES CONTOURS
Minip=0.5;"Minimum planet mass (Msat)";
Maxip=3.0;"Maximum planet mass (Msat)";

Minim=1.0;"Minimum moon mass (Mencel)";
Maxim=12.0;"Maximum moon mass (Mencel)";

#VECTORS FOR PLANET AND MOON MASSES
N=100.0
Mps=np.linspace(Minip,Maxip,N)
Mms=np.linspace(Minim,Maxim,N)
MMS,MPS=np.meshgrid(Mms,Mps)

#ROUTINE'S ARGUMENTS
Rp=MPS**(0.156)
Ms=Msu*Msun
am=amu*Rp
Mp=MPS*Msat
Mm=MMS*Mencel
P=tsu/YEAR

#FINALLY WE GET THE ROUND-TRIP TIMESCALE

tsync=syncTime(tini,numiter,P,Ms,Mp,Mm,am,A)
tdecay=decayTime(P,Ms,Mp,Mm,am,tsync,A)
trt=(tsync+tdecay)/Gyear



print (trt.max(),trt.min())

fig=plt.figure(figsize=(7.5,7.5))

inf=trt.min()+0.00001
sup=trt.max()+0.00001

levels=np.linspace(inf,sup,100)
ct=plt.contourf(MMS,MPS,trt,levels=levels,coloring='heatmap',autocontour=False,erase=True,title="right")
cbar=fig.colorbar(ct,format="%.2f")
cbar.set_label('Round-trip times')

#cbar.set_label('Round-trip times', rotation=270)
#colorbar=(title='Round-trip times',titleside='right',titlefont=Font(size=14,Family='Arial, sans-serif'))

levels=np.linspace(inf,sup,10.0)
ct=plt.contour(MMS,MPS,trt,levels=levels,colors='k',labels='inline')
plt.clabel(ct,inline=1,fmt="%.2f Gyr")
plt.ylabel("$M_{\mathrm{p}}$ / $M_{\mathrm{Sat}}$",fontsize=15)
plt.xlabel("$M_{\mathrm{m}}$ / $M_{\mathrm{Enceladus}}$",fontsize=15)
plt.xticks(fontsize=13.0)
plt.yticks(fontsize=13.0)

#THE PLANET SEMIMAJOR AXIS ap
np=2*pi/P
ap=((G*Ms)/np**(2.0))**(1.0/3.0)

print ((ap*6371000)/1.5e11)

plt.title(r"$a_{\mathrm{m,ini}}=%.2f\,R_{\mathrm{p}}$ , $P_{\mathrm{orb}}=%.2f$ d , $\alpha=%.2f$"%(amu,tsu,A))
plt.savefig("Mass-contours.png")



###########CONTOURS FOR THE PLANETARY AND MOON POSITIONS##########
##################################################################


#FIXED PARAMETERS
Msu=1.0;"Stellar mass (Msun)";
Mpu=1.0;"Planet mass (Msat)";
Mmu=10.0;"Moon mass (Mencel)";
A=0.18;"Alpha of the planet";

#ROUTINE'S ARGUMENTS
Rp=Mpu**(0.156)
Ms=Msu*Msun
Mp=Mpu*Msat
Mm=Mmu*Mencel

#CALCULUS OF K2/Q

eps2=(Oi/(G*Mp/Rp**3.0)**(0.5))**2.0
B=kpp*A**bp*(Mp*95.159)**delta;"Beta of the planet";
Y=(A**3.0*(1-B))/(B*(1-A**3.0));"Gamma of the planet";
T1=100*pi/63
T2=A**5.0/(1-A**5.0)
T3=1+((1-Y)/Y)*A**3.0
T4=(1+(2.5*(1-Y)/Y)*A**3.0)**(-2.0)
KQ=T1*eps2*T2*T3*T4

print (KQ)

#TO KNOW THE ROCHE RADIUS
DensP=(3.0*Mp)/(4*pi*Rp**(3.0))
B1=(DensP/DensPart)**(1.0/3.0)
ar=2.0*B1*Rp

print ("Roche Radius= %.3f Rp"%(ar))

nsyncini=Oi
asyncini=((G*Mp)/nsyncini**(2.0))**(1.0/3.0)
asyncini=(asyncini/Rp)

print ("Initial synchronous radius= %.3f Rp"%(asyncini))
#Cuando ar sea 2.5 y asyncini 2.7 se debe escoger 2.7
amin=max(ar,asyncini)

#FOR THE POSITIONS CONTOURS
tinip=20.0;"Minimum orbital period of the planet (days)";
tmaxp=100.0;"Maximum orbital period of the planet (days)";

ainim=1.1*amin;"Minimum moon position (Rsat)";
amaxim=2*amin;"Maximum moon position (Rsat)";

#VECTORS FOR PLANET ORBITAL PERIOD AND MOON POSITIONS
tps=np.linspace(tinip,tmaxp,100)
ams=np.linspace(ainim,amaxim,100)
AMS,TPS=np.meshgrid(ams,tps)

P=TPS/YEAR
am=AMS*Rp

#FINALLY WE GET THE ROUND-TRIP TIMESCALE
tsync=syncTime(tini,numiter,P,Ms,Mp,Mm,am,A)
tdecay=decayTime(P,Ms,Mp,Mm,am,tsync,A)
trt=(tsync+tdecay)/Gyear

fig=plt.figure(figsize=(7.5,7.5))

inf=trt.min()+0.00001
sup=trt.max()+0.00001

levels=np.linspace(inf,sup,100)
ct=plt.contourf(AMS,TPS,trt,levels=levels,coloring='heatmap',autocontour=False,erase=True)
cbar=fig.colorbar(ct,format="%.2f")
cbar.set_label('Round-trip times')
levels=np.linspace(inf,sup,10.0)
ct=plt.contour(AMS,TPS,trt,levels=levels,colors='k',labels='inline')
plt.clabel(ct,inline=1,fmt="%.2f Gyr")
plt.ylabel("$P_{\mathrm{orb}}\,[\mathrm{days}]$",fontsize=15)
plt.xlabel("$a_{\mathrm{m,ini}}$ / $R_{\mathrm{p}}$",fontsize=15)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.title(r"$M_{\mathrm{p}}=%.2f\,M_{\mathrm{Sat}} $ , $M_{\mathrm{m}}=%.1f\,M_{\mathrm{Enceladus}}$, $\alpha=%.2f$"%(Mpu,Mmu,A))

plt.savefig("Position-contours.png")




#############THE FITS FOR THE DISSIPATION FUNCTION (EMPIRICAL)#############
###########################################################################

import numpy as np
from math import*
from matplotlib import pyplot as plt
from sys import exit
tidal=np.loadtxt("Tidalscale.csv",skiprows=1,delimiter=";",dtype='|S9')
tidal1=np.loadtxt("Tidalscale-kip.csv",skiprows=1,delimiter=";",dtype='|S9')

radius=tidal[:,1]
radius=radius.astype(np.float)

ratiok2Q=tidal[:,2]
ratiok2Q=ratiok2Q.astype(np.float)

data=[]
for i in range(len(ratiok2Q)):
    data+=[[radius[i],ratiok2Q[i]]]
np.savetxt("data.dat",data)

setvalues=np.loadtxt("data.dat")

radius1=tidal1[:,1]
radius1=radius1.astype(np.float)

ratiok2Q1=tidal1[:,2]
ratiok2Q1=ratiok2Q1.astype(np.float)

data1=[]
for i in range(len(ratiok2Q1)):
    data1+=[[radius1[i],ratiok2Q1[i]]]
np.savetxt("data1.dat",data1)

setvalues1=np.loadtxt("data1.dat")

#DESDE AQUI COMIENZA EL AJUSTE LINEAL PARA DATA

RP=setvalues[:,0]
TD=setvalues[:,1]

LogRP=np.log10(RP)
LogTD=np.log10(TD)

z=np.polyfit(LogRP,LogTD,1)
zz=np.poly1d(z)
ys=zz(LogRP)

m=zz[1]
b=zz[0]

#DESDE AQUI COMIENZA EL AJUSTE LINEAL PARA DATA1

RP1=setvalues1[:,0]
TD1=setvalues1[:,1]

LogRP1=np.log10(RP1)
LogTD1=np.log10(TD1)

w=np.polyfit(LogRP1,LogTD1,1)
ww=np.poly1d(w)
zs=ww(LogRP1)

x=ww[1]
xx=ww[0]

'''
bup=-2.37
bdown=-3.67
Upper=zz[1]*LogRP+bup
Lower=zz[1]*LogRP+bdown
'''

print(10**(m*np.log10(9.14)+b))
print(b)
#print (zz)

plt.figure(figsize=(8.0,8.0))
plt.plot(LogRP,LogTD,'go',label='Observed values')
plt.plot(LogRP1,LogTD1,'bo')
plt.plot(LogRP,ys,'k-',label='Linear fitting')
plt.plot(LogRP1,zs,'b-',label='Without Saturn')
#plt.plot(LogRP,Upper,'r--',label='Upper bound')
#plt.plot(LogRP,Lower,'b--',label='Lower bound')
plt.title("Empirical relation for $k_{2}/Q$")
plt.legend(loc='best',fontsize=14,frameon=False)
plt.xlabel(r"$\mathrm{Log}(R_{\mathrm{p}}/R_{\oplus})$",weight='medium',fontsize=16)
plt.ylabel("$\mathrm{Log}(k_{2}/Q)$",weight='medium',fontsize=15)
plt.xlim((min(LogRP)-0.1,max(LogRP)+0.1))
plt.savefig("plotvaluesk2Q.png")


##THIS IS FOR THE CONFIRMED PLANETS AND KOIS, SINCE THERE IS NOT Mp BUT Rp###
#############################################################################


def g(t,P,Ms,Rp,Mm,am,A):
    np=(2*pi)/P
    ap=((G*Ms)/np**(2.0))**(1.0/3.0)
    Mp=(Rp)**(1.0/0.156)*Msat
    eps2=(Oi/(G*Mp/Rp**3.0)**(0.5))**2.0
    B=kpp*A**bp*(Mp*95.159)**delta;"Beta of the planet";
    Y=(A**3.0*(1-B))/(B*(1-A**3.0));"Gamma of the planet";
    T1=100*pi/63
    T2=A**5.0/(1-A**5.0)
    T3=1+((1-Y)/Y)*A**3.0
    T4=(1+(2.5*(1-Y)/Y)*A**3.0)**(-2.0)
    KQ=T1*eps2*T2*T3*T4
    nmi=((G*(Mp))/am**(3.0))**(0.5)
    A1=(nmi**(-13.0/3.0)+(19.5*KQ*Mm*(Rp)**(5.0)*t)/(G**(5.0/3.0)*(Mp)**(8.0/3.0)))
    A2=(Oi-(A1)**(-3.0/13.0))
    A3=A2*((2*GR*(Mp))/(3*KQ*G*Rp**(3.0)))
    A4=(Ms**(2.0)/ap**(6.0)+(Mm)**(2.0)/am**(6.0))
    Gfun=A3*A4**(-1.0)
    return Gfun

#THIS IS THE FIXED-POINT FUNCTION FOR DOING THE "k" ITERATIONS TO
#REFINE THE SEED'S VALUE (tini).

def syncTime(t,k,P,Ms,Rp,Mm,am,A):
    for i in np.arange(k):
        t=g(t,P,Ms,Rp,Mm,am,A)
    return t

#THIS IS GOING TO BE THE DATA FOR DOING THE CALCULUS OF nsync AND async

def decayTime(P,Ms,Rp,Mm,am,tsync,A):
    np=(2*pi)/P
    ap=((G*Ms)/np**(2.0))**(1.0/3.0)
    Mp=(Rp)**(1.0/0.156)*Msat
    eps2=(Oi/(G*Mp/Rp**3.0)**(0.5))**2.0
    B=kpp*A**bp*(Mp*95.159)**delta;"Beta of the planet";
    Y=(A**3.0*(1-B))/(B*(1-A**3.0));"Gamma of the planet";
    T1=100*pi/63
    T2=A**5.0/(1-A**5.0)
    T3=1+((1-Y)/Y)*A**3.0
    T4=(1+(2.5*(1-Y)/Y)*A**3.0)**(-2.0)
    KQ=T1*eps2*T2*T3*T4
    nmi=((G*(Mp))/am**(3.0))**(0.5)
    NMI=nmi**(-13.0/3.0)
    
    #CALCULUS OF THE ROCHE RADIUS
    DensP=(3.0*(Mp))/(4*pi*Rp**(3.0))
    B1=(DensP/DensPart)**(1.0/3.0)
    ar=2.0*B1*Rp
    #print ("Roche Radius= %.2f Rsat"%(ar))
    
    NUM=(19.5*KQ*(Mm)*Rp**(5.0)*tsync)
    DEN=(G**(5.0/3.0)*(Mp)**(8.0/3.0))
    nsync=(NMI+NUM/DEN)**(-3.0/13.0)
    async=((G*(Mp))/nsync**(2.0))**(1.0/3.0)

    #NEXT, WE CALCULATE THE DECAY TIME FROM THE SYNCHRONOUS RADIUS TO THE
    #ROCHE RADIUS
    TER1=(19.5*KQ*(Mm)*(Rp)**(5.0))
    TER2=(G**(5.0/3.0)*(Mp)**(8.0/3.0))
    TER3=TER2/TER1
    TER4=(-((G*(Mp))/ar**(3.0))**(-13.0/6.0)+nsync**(-13.0/3.0))
    Td=TER3*TER4
    return Td


##################AT FIRST FOR THE CONFIRMED PLANETS########################
############################################################################

#HERE WE OBTAIN THE VALUES OF Mp AND Porb FOR THE CONFIRMED PLANETS BY KEPLER
filtered=np.loadtxt("Exoplanets_filtered1.csv",skiprows=1,delimiter=";",dtype='|S9')

Rpp=filtered[:,2]
Rpp=Rpp.astype(np.float)

PP=filtered[:,4]
PP=PP.astype(np.float)

#FIXED PARAMETERS
Msu=1.0;"Stellar mass (Msun)";
amu=4.0;"Initial moon position (Rp)";
mmu=10.0;"Mass of the moon (Mencel)";

#ROUTINE PARAMETERS

Rp=Rpp
Ms=Msu*Msun
am=amu*Rp
Mm=mmu*Mencel
P=PP/YEAR

#FINALLY WE GET THE ROUND-TRIP TIMESCALE

Amin=0.11
Amax=0.23
Ahere=Amin+(Amax-Amin)*np.random.random(100)
tsync=np.array([syncTime(tini,numiter,P,Ms,Rp,Mm,am,A) for A in Ahere])
tdecay=np.array([decayTime(P,Ms,Rp,Mm,am,t,A) for A,t in zip(Ahere,tsync)])
trt=(tsync+tdecay)/Gyear

trtm=np.array([np.mean(trt[:,i]) for i in range(trt.shape[1])])
trts=np.array([np.std(trt[:,i]) for i in range(trt.shape[1])])

trtmin=np.array([min(trt[:,i]) for i in range(trt.shape[1])])
trtmax=np.array([max(trt[:,i]) for i in range(trt.shape[1])])

trt=trtm


round_trip_confirmedp=[]
for i in range(len(PP)):
    round_trip_confirmedp+=[[trt[i],trts[i]]]
np.savetxt("round_trip_confirmedp.csv",round_trip_confirmedp)

print (round_trip_confirmedp)

print (trt.max(),trt.min(),len(PP))

#CURVAS DE RADIO DE HILL EN RADIOS PLANETARIOS COMO FUNCION DEL PERIODO ORBITAL DEL PLANETA.
#TERMINOS FIJOS 
DAY=86400.
RSAT=58232000.
MSAT=5.6836E26
bet=(1.0/3.0)-0.156
f=0.48

#CALCULO DE LA CONSTANTE
NUM=(DAY**2.0)*6.67384e-11
DEN=12.0*pi**2.0
TERM=(NUM/DEN)**(1.0/3.0)
TERM1=MSAT**(1.0/3.0)/RSAT
Cons=TERM*TERM1

#ESTE ES EL PERIODO EN DIAS
Peri=np.linspace(0.1,100.0,100.0)

plt.figure(figsize=(7.5,7.5))

hard=0.3

MP_old=np.zeros(len(Peri))
for RH in np.arange(5.0,100.0,5.0):
    MP=(1.0/Cons)**(1.0/bet)*(RH/f)**(1.0/bet)*Peri**(-2.0/(3.0*bet))
    color=cm.jet(1.0-(RH/90.0)**0.6)
    plt.fill_between(Peri,MP_old,MP,facecolor=color,alpha=0.5)
    MP_old=MP
    

for i in np.arange(len(trt)):
    if trt[i]<0.05:
        deg=4
        col='b'
    elif trt[i]<0.1:
        deg=6
        col='b'
    elif trt[i]<0.3:
        deg=9
        col='b'
    elif trt[i]<0.5:
        deg=11
        col='b'
    elif trt[i]<1.0:
        deg=14
        col='b'
    elif trt[i]<1.5:
        deg=17
        col='b'
    elif trt[i]<2.0:
        deg=21
        col='b'
    elif trt[i]<2.5:
        deg=25
        col='b'
    elif trt[i]<3.0:
        deg=27
        col='b'
    else:
        deg=32
        col='r'
    plt.plot(PP[i],Rpp[i],'o',color=col,ms=deg,markeredgecolor='none',alpha=hard,zorder=100)

    deg=trt[i]-trts[i]
    if deg>3.0:
        deg=32
        color='r'
    else:
        color='b'
        deg=deg*2.4
    plt.plot(PP[i],Rpp[i],'o',markerfacecolor='none',ms=deg,markeredgecolor=color
             ,zorder=100)
    
    deg=trt[i]+trts[i]
    if deg>3.0:
        deg=32
        color='r'
    else:
        color='b'
        deg=deg
    plt.plot(PP[i],Rpp[i],'o',markerfacecolor='none',ms=deg,markeredgecolor=color,
             zorder=100)

dob=20.0
Trt1=plt.scatter([],[],marker='o',color='k',alpha=0.6,edgecolors='none',s=dob)
Trt2=plt.scatter([],[],marker='o',color='k',alpha=0.6,edgecolors='none',s=4*dob)
Trt3=plt.scatter([],[],marker='o',color='k',alpha=0.6,edgecolors='none',s=13*dob)
Trt4=plt.scatter([],[],marker='o',color='k',alpha=0.6,edgecolors='none',s=30*dob)
Trt5=plt.scatter([],[],marker='o',color='k',alpha=0.6,edgecolors='none',s=50*dob)
labels=["","","","",""]
legend=plt.legend([Trt1,Trt2,Trt3,Trt4,Trt5],labels,ncol=5,frameon=True,fontsize=10.0,handlelength=1.0,loc='upper right',borderpad=1.25,handletextpad=0.5,title=r'$0.0\,<\,t_{\mathrm{round-trip}}\,[\mathrm{Gyr}]<\,3.0$',scatterpoints=1)
frame=legend.get_frame()
frame.set_facecolor('0.95')
legend.set_zorder(105)
plt.ylabel("$R_{\mathrm{p}}$ / $R_{\mathrm{Sat}}$",fontsize=15)
plt.xlabel("$P_{\mathrm{orb}}\,[\mathrm{days}]$",fontsize=15)
plt.xticks(fontsize=13.0)
plt.yticks(fontsize=13.0)
plt.ylim((0.6,1.4))
plt.title(r"$a_{\mathrm{m,ini}}=%.1f\,R_{\mathrm{p}}$, $M_{\mathrm{m}}=%.1f\,M_{\mathrm{Enceladus}}$"%(amu,mmu))
plt.grid()
plt.savefig("populationwarmplanets-mean.png")


########################SECONDLY FOR THE KOIS###############################
############################################################################

#HERE WE OBTAIN THE VALUES OF Mp AND Porb FOR THE CONFIRMED PLANETS BY KEPLER
filtered=np.loadtxt("cumulative.csv",skiprows=1,delimiter=";",dtype='|S9')

Rpp=filtered[:,6]
Rpp=Rpp.astype(np.float)

PP=filtered[:,4]
PP=PP.astype(np.float)

#FIXED PARAMETERS
Msu=1.0;"Stellar mass (Msun)";
amu=4.0;"Initial moon position (Rp)";
mmu=10.0;"Mass of the moon (Mencel)";

#ROUTINE PARAMETERS

Rp=Rpp
Ms=Msu*Msun
am=amu*Rp
Mm=mmu*Mencel
P=PP/YEAR

#FINALLY WE GET THE ROUND-TRIP TIMESCALE

Amin=0.11
Amax=0.23
Ahere=Amin+(Amax-Amin)*np.random.random(100)
tsync=np.array([syncTime(tini,numiter,P,Ms,Rp,Mm,am,A) for A in Ahere])
tdecay=np.array([decayTime(P,Ms,Rp,Mm,am,t,A) for A,t in zip(Ahere,tsync)])
trt=(tsync+tdecay)/Gyear

trtm=np.array([np.mean(trt[:,i]) for i in range(trt.shape[1])])
trts=np.array([np.std(trt[:,i]) for i in range(trt.shape[1])])

trtmin=np.array([min(trt[:,i]) for i in range(trt.shape[1])])
trtmax=np.array([max(trt[:,i]) for i in range(trt.shape[1])])

trt=trtm

round_trip_kois=[]
for i in range(len(PP)):
    round_trip_kois+=[[trt[i],trts[i]]]
np.savetxt("round_trip_kois.csv",round_trip_kois)

print (trt.max(),trt.min(),len(PP))

#CURVAS DE RADIO DE HILL EN RADIOS PLANETARIOS COMO FUNCION DEL PERIODO ORBITAL DEL PLANETA.
#TERMINOS FIJOS 
DAY=86400.
RSAT=58232000.
MSAT=5.6836E26
bet=(1.0/3.0)-0.156
f=0.48

#CALCULO DE LA CONSTANTE
NUM=(DAY**2.0)*6.67384e-11
DEN=12.0*pi**2.0
TERM=(NUM/DEN)**(1.0/3.0)
TERM1=MSAT**(1.0/3.0)/RSAT
Cons=TERM*TERM1

#ESTE ES EL PERIODO EN DIAS
Peri=np.linspace(0.1,100.0,100.0)

#ESTE ES EL CALCULO DE LA MASA PLANETARIA EN MASAS DE SATURNO

plt.figure(figsize=(7.5,7.5))

hard=0.3

MP_old=np.zeros(len(Peri))
for RH in np.arange(5.0,100.0,5.0):
    MP=(1.0/Cons)**(1.0/bet)*(RH/f)**(1.0/bet)*Peri**(-2.0/(3.0*bet))
    color=cm.jet(1.0-(RH/90.0)**0.6)
    plt.fill_between(Peri,MP_old,MP,facecolor=color,alpha=0.5)
    MP_old=MP
    

for i in np.arange(len(trt)):
    if trt[i]<0.05:
        deg=4
        col='b'
    elif trt[i]<0.1:
        deg=6
        col='b'
    elif trt[i]<0.3:
        deg=9
        col='b'
    elif trt[i]<0.5:
        deg=11
        col='b'
    elif trt[i]<1.0:
        deg=14
        col='b'
    elif trt[i]<1.5:
        deg=17
        col='b'
    elif trt[i]<2.0:
        deg=21
        col='b'
    elif trt[i]<2.5:
        deg=25
        col='b'
    elif trt[i]<3.0:
        deg=27
        col='b'
    else:
        deg=32
        col='r'
    plt.plot(PP[i],Rpp[i],'o',color=col,ms=deg,markeredgecolor='none',alpha=hard,zorder=100)

    deg=trt[i]-trts[i]
    if deg>3.0:
        deg=32
        color='r'
    else:
        color='b'
        deg=deg*1.0
    plt.plot(PP[i],Rpp[i],'o',markerfacecolor='none',ms=deg,markeredgecolor=color
             ,zorder=100)
    
    deg=trt[i]+trts[i]
    if deg>3.0:
        deg=32
        color='r'
    else:
        color='b'
        deg=deg
    plt.plot(PP[i],Rpp[i],'o',markerfacecolor='none',ms=deg,markeredgecolor=color,
             zorder=100)

dob=20.0
Trt1=plt.scatter([],[],marker='o',color='k',alpha=0.6,edgecolors='none',s=dob)
Trt2=plt.scatter([],[],marker='o',color='k',alpha=0.6,edgecolors='none',s=4*dob)
Trt3=plt.scatter([],[],marker='o',color='k',alpha=0.6,edgecolors='none',s=13*dob)
Trt4=plt.scatter([],[],marker='o',color='k',alpha=0.6,edgecolors='none',s=30*dob)
Trt5=plt.scatter([],[],marker='o',color='k',alpha=0.6,edgecolors='none',s=50*dob)
labels=["","","","",""]
legend=plt.legend([Trt1,Trt2,Trt3,Trt4,Trt5],labels,ncol=5,frameon=True,fontsize=10.0,handlelength=1.0,loc='upper right',borderpad=1.25,handletextpad=0.5,title=r'$0.0\,<\,t_{\mathrm{round-trip}}\,[\mathrm{Gyr}]<\,3.0$',scatterpoints=1)
frame=legend.get_frame()
frame.set_facecolor('0.95')
legend.set_zorder(105)
plt.ylabel("$R_{\mathrm{p}}$ / $R_{\mathrm{Sat}}$",fontsize=15)
plt.xlabel("$P_{\mathrm{orb}}\,[\mathrm{days}]$",fontsize=15)
plt.xticks(fontsize=13.0)
plt.yticks(fontsize=13.0)
plt.ylim((0.8,1.4))
plt.title(r"$a_{\mathrm{m,ini}}=%.1f\,R_{\mathrm{p}}$, $M_{\mathrm{m}}=%.1f\,M_{\mathrm{Enceladus}}$"%(amu,mmu))
plt.grid()
plt.savefig("koiscandidates.png")



