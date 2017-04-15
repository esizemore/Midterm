#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 11:52:42 2017

@author: elizabethsizemore
"""

import numpy as np
import matplotlib.pyplot as plt
from math import *


#Constants
ro_Pb=11.3*1000 #kg/m^3 #density of Pb
ro_Ba= 3.6*1000 #kg/m^3 #density of Ba
ro_Np= 20.2 *1000 #kg/m^3 #density of Np
r=0.04 #radius in m
theta=48*np.pi/180 #48 degrees converted to radians
x0=0 #initial position in x direction (m)
y0=5 #initial position in y direction (m)
v=120.0 #initial velocity m/s
vx=v*np.cos(theta) #x component of velocity
vy=v*np.sin(theta) #y component of velocity
ro_air=1.26 #kg/m^3
C= 0.47 #coefficient of drag
g=9.81 



#PROGRAM FOR PB CANNONBALL

m_Pb=(4/3)*np.pi*(r**3)*ro_Pb #gives mass of spherical Pb cannonball


#split two second order eqns into four first order eqns
def fPb(r_Pb, t): #define function wrt position and time
    Pbx=r_Pb[0] #position
    Pby=r_Pb[1] #position
    Pbvx=r_Pb[2] #velocity
    Pbvy=r_Pb[3] #velocity
    Pbax= (-1*np.pi*(r**2)*ro_air* C * Pbvx*v)/(2*m_Pb) #acceleration
    Pbay=((-1)*g)- ((np.pi *(r**2))*ro_air* C * Pbvy*v)/(2*m_Pb) #acceleration
    return np.array([Pbvx, Pbvy, Pbax, Pbay], float)
    
    
    
a = 0.0
b= 1000.0
N=100000
h=(b-a)/N


tpoints=np.arange(a,b,h)
Pb_xpoints =[]
Pb_ypoints =[]

r_Pb= np.array ([x0,y0,vx,vy],float)


for t in tpoints:
    Pb_xpoints.append(r_Pb[0])
    Pb_ypoints.append(r_Pb[1])
    k1=h*fPb(r_Pb, t)
    k2= h*fPb(r_Pb+ 0.5*k1, t+0.5*h)
    k3=h*fPb(r_Pb + 0.5*k2, t+0.5*h)
    k4=h*fPb(r_Pb+k3, t+h)
    r_Pb+= (k1+2*k2+2*k3+k4)/6
    if r_Pb[1] < 0:
        break
    else:
        continue
    
    

    
plt.plot(Pb_xpoints, Pb_ypoints, label='Pb')
plt.title('Projectile Motion of Cannonball')
plt.xlabel('Distance (m)')
plt.ylabel('Height (m)')
plt.show()
    










#PROGRAM FOR BA CANNONBALL

m_Ba=(4/3)*np.pi*(r**3)*ro_Ba #gives mass of spherical Pb cannonball




#split two second order eqns into four first order eqns
def fBa(r_Ba, t): #define function wrt position and time
    Bax=r_Ba[0] #position
    Bay=r_Ba[1] #position
    Bavx=r_Ba[2] #velocity
    Bavy=r_Ba[3] #velocity
    Baax= (-1*np.pi*(r**2)*ro_air* C * Bavx*v)/(2*m_Ba) #acceleration
    Baay=((-1)*g)- ((np.pi *(r**2))*ro_air* C * Bavy*v)/(2*m_Ba) #acceleration
    return np.array([Bavx, Bavy, Baax, Baay], float)
    
    
    
a = 0.0
b= 1000.0
N=100000
h=(b-a)/N


tpoints=np.arange(a,b,h)
Ba_xpoints =[]
Ba_ypoints =[]

r_Ba= np.array ([x0,y0,vx,vy],float)


for t in tpoints:
    Ba_xpoints.append(r_Ba[0])
    Ba_ypoints.append(r_Ba[1])
    k1=h*fBa(r_Ba, t)
    k2= h*fBa(r_Ba+ 0.5*k1, t+0.5*h)
    k3=h*fBa(r_Ba + 0.5*k2, t+0.5*h)
    k4=h*fBa(r_Pb+k3, t+h)
    r_Ba+= (k1+2*k2+2*k3+k4)/6
    if r_Ba[1] < 0:
        break
    else:
        continue
    
    

    
plt.plot(Ba_xpoints, Ba_ypoints, label='Ba')
plt.title('Projectile Motion of Cannonball')
plt.xlabel('Distance (m)')
plt.ylabel('Height (m)')
plt.show()
    
    
    
    
    
    
    


#PROGRAM FOR NP CANNONBALL 
    
    
m_Np=(4/3)*np.pi*(r**3)*ro_Np #gives mass of spherical Pb cannonball




#split two second order eqns into four first order eqns
def fNp(r_Np, t): #define function wrt position and time
    Npx=r_Np[0] #position
    Npy=r_Np[1] #position
    Npvx=r_Np[2] #velocity
    Npvy=r_Np[3] #velocity
    Npax= (-1*np.pi*(r**2)*ro_air* C * Npvx*v)/(2*m_Np) #acceleration
    Npay=((-1)*g)- ((np.pi *(r**2))*ro_air* C * Npvy*v)/(2*m_Np) #acceleration
    return np.array([Npvx, Npvy, Npax, Npay], float)
    
    
    
a = 0.0
b= 1000.0
N=100000
h=(b-a)/N


tpoints=np.arange(a,b,h)
Np_xpoints =[]
Np_ypoints =[]

r_Np= np.array ([x0,y0,vx,vy],float)


for t in tpoints:
    Np_xpoints.append(r_Np[0])
    Np_ypoints.append(r_Np[1])
    k1=h*fNp(r_Np, t)
    k2= h*fNp(r_Np+ 0.5*k1, t+0.5*h)
    k3=h*fNp(r_Np + 0.5*k2, t+0.5*h)
    k4=h*fNp(r_Np+k3, t+h)
    r_Np+= (k1+2*k2+2*k3+k4)/6
    if r_Np[1] < 0:
        break
    else:
        continue
    
    

    
plt.plot(Np_xpoints, Np_ypoints, label='Np')
plt.title('Projectile Motion of Cannonball')
plt.xlabel('Distance (m)')
plt.ylabel('Height (m)')
plt.show()
    
    
    
print ('Maximum Distance with Pb cannonball', Pb_xpoints[-1])
print ('Maximum Distance with Ba cannonball', Ba_xpoints[-1])
print ('Maximum Distance with Np cannonball', Np_xpoints[-1])

    






#OUTPUT OF FILE
#>> runfile('/Users/elizabethsizemore/Documents/Computational_Physics/Sizemore_Midterm.py', wdir='/Users/elizabethsizemore/Documents/Computational_Physics')
#Maximum Distance with Pb cannonball 826.65418088
#Maximum Distance with Ba cannonball 395.878718415
#Maximum Distance with Np cannonball 1031.98073271
#>>> 


#The Np cannonball travels the furthest horizontal distance. This 
#makes sense because Np had a greater density






