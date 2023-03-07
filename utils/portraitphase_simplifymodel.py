#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 17:47:20 2023

@author: loms
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from scipy import integrate
# import scienceplots

# plt.style.use(['science','no-latex'])
plt.close('all')

# Set up

umax1=1.9872  # maximum growth rates of P1 [d**{-1}]
umax2=2.7648  # maximum growth rates of P2 [d**{-1}]
gmax1=0.3     # maximum grazing rates of Z on P1 [d**{-1}]
gmax2=0.6     # maximum grazing rates of Z on P2 [d**{-1}]
mP=0		  # P2 mortality rate (default 0 ie no P2 sinking) [d**{-1}]
mZ=0.05 	  # Z2 quadratic mortality rate [mmolC**{-1} m**{3} d**{-1}]
eZ=0.1        # zoo excretion rate (Z1 and Z2) [d**{-1}]
gamma1=0.6    # conversion factor from P1 to Z [/]
gamma2=0.6    # conversion factor from P2 to Z [/]
epsilon=0.25  # fraction of Z excretion that is available as regenerated PO4 [/]
Psupply = 0.01

def growth_model_simplify(u,t=0):

    return ([ u[3]*umax1*u[0]-u[0]*gmax1*u[2], u[3]*umax2*u[1]-u[1]*gmax2*u[2], gamma1*u[0]*gmax1*u[2]+gamma2*u[1]*gmax2*u[2]-eZ*((1-gamma1)*u[0]*gmax1*u[2]+(1-gamma2)*u[1]*gmax2*u[2])-mZ*u[2]**2, 
           Psupply+epsilon*eZ*((1-gamma1)*u[0]*gmax1*u[2]+(1-gamma2)*u[1]*gmax2*u[2])-(u[3]*umax1*u[0]+u[3]*umax2*u[1]) ])

""" Equilibre du modèle calculés analytiquement """

# 1/ P1 = 0 
print('Equilibrium values for P1 =0 :')
PO4eq1 = Psupply*(gmax2*(gamma2-eZ*(1-gamma2))*(epsilon*eZ*(1-gamma2)+1))/mZ
print('PO4 =', PO4eq1)
Zeq1 = (umax2/gmax2)*PO4eq1
print('Z =', Zeq1)
P2eq = mZ*Zeq1/(gamma2*gmax2-eZ*(1-gamma2*gmax2))
print('P2 =',P2eq)

EQ2 = ([0,P2eq, Zeq1, PO4eq1])

# 2/ P2 = 0
print('Equilibrium values for P2 =0 :')
PO4eq2 = Psupply*(gmax1*(gamma1-eZ*(1-gamma1))*(epsilon*eZ*(1-gamma1)+1))/mZ
PO4eq2 = np.array(PO4eq2)
print('PO4 =',PO4eq2)
Zeq2 = (umax1/gmax1)*PO4eq2
Zeq2 = np.array(Zeq2)
print('Z =',Zeq2)
P1eq = mZ*Zeq2/(gamma1*gmax1-eZ*(1-gamma1*gmax1))
P1eq = np.array(P1eq)
print('P1 =',P1eq)
      
EQ1 = ([P1eq,0, Zeq2,PO4eq2])          

"""  Trajectoires """
valuesP = np.linspace(0.1,5,10)
vcolors=plt.cm.autumn_r(np.linspace(0.1,1,len(valuesP)))

#P1 = 0
plt.figure(1)
ax = plt.axes(projection ='3d')
t = np.linspace(0,150,1000)

for i, col1 in zip(valuesP, vcolors):
            u0 = [E*i for E in EQ2]
            u = integrate.odeint(growth_model_simplify,u0,t)
            ax.plot(u[:,1], u[:,2], u[:,3], lw=0.5*i, color=col1)
            ax.view_init(23, -151)
             
# P2 = 0
plt.figure(2)
ax = plt.axes(projection ='3d')
t = np.linspace(0,150,1000)

for i, col1 in zip(valuesP, vcolors):
            u0 = [E*i for E in EQ1]
            u = integrate.odeint(growth_model_simplify,u0,t)
            ax.plot(u[:,0], u[:,2], u[:,3], lw=0.5*i, color=col1)
            ax.view_init(23, -151)          
            
            