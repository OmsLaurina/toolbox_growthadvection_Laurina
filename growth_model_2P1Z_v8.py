#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 11:23:36 2023

@author: loms
"""
"""
Implementation of the Growth Model 
                iP1Z
with growth rate and grazing rates = monod equation

i1 = 1 phytoplancton
i2 = 2 phytoplankton
                
Created on Tue Dec 20 12:39:07 2022

@author: loms
"""
import matplotlib.pyplot as plt
import numpy as npy
import numpy as np

plt.close('all')

# timestep
dt = 0.2
# the time at which the simulation ends (days)
end_time = 600
# creates a time vector from 0 to end_time, seperated by a timestep
time = npy.arange(0,end_time,dt)

"""" definition of model variables"""
##Creating empty vectors
# growth rates (d^-1)
u_1=[]
u_2=[]			
# grazing rates (d^-1)
g_1=[]
g_2=[]				
#Primary production
PP_1=[]
PP_2=[]
#Grazing
G1=[]
G2=[]
# biomass of Phytoplankton and Zooplankton (mmolC/m³)
P_1=[] 			
P_2=[] 						
Z = []
#Phosphate concentration
PO4 = []
Psupply=[]
#exc and regeneration
exc_Z=[]
PO4_reg=[]
p_fec=[]
#death
death_Z=[]
#export
Export=[]

# initial conditions at time=0	
P0_1 = 0.6
P0_2 = 0.1
Z0 = 0.6
PO40 = 0.5

"""" definition of model parameters"""

## Monod parameters
# growth rates (d^-1)
umax_1 = 1.9872
umax_2 = 2.7648 
# grazing rates (d^-1)
gmax_1 = 1.4226
gmax_2= 1.4226
# half saturation constants
kP_1 = 1
kP_2 = 2
kZ = 5

## Other model parameters
# Z exc fraction (d^-1)
eZ = 0.1
# fraction of Z exc that is available as regenerated PO4
epsilon = 0.75
# quadratic mortality of Z (mmolC^{-1} m^{-3} d^{-1})
mZ = 0.005

##Apport externe de nutriment 
Psupply_moy = 0.02

#Psupply constant
Psupply =[Psupply_moy] * len(time)

#Psupply variable

#Sinusoïdalement
# b = 0.05
# w = 1
# Psupply_sin = []
# for i in range(len(time)):
#     Psupply_sin.append(Psupply_moy+ b * np.sin(w*i))
#     Psupply = Psupply_sin 

#Pulsé
# T = 4
# b = (Psupply_moy*0.25*10)/(1-np.exp(-0.25*T))-Psupply_moy
# Psupply_sin = []
# for i in range(len(time)):
#     Psupply_sin.append(Psupply_moy+b*np.sin(w*i))
#     Psupply = Psupply_sin 



def ga_model_2P1Z_v8():
    
    ## -------------- Parameters
    global dt, end_time, time, umax_1, umax_2, gmax_1, gmax_2, kP_1, kP_2, KZ, eZ, epsilon, mZ, u_1, u_2, P_1, P_2, Z, u0_1, u0_2, P0_1, P0_2, Z0, Psupply
    
    # Initial conditions at time=0			
    P_1.append(P0_1)
    P_2.append(P0_2)
    Z.append(Z0)
    PO4.append(PO40)
    
    
    nb_time=len(time)
    
    for t in range(0,nb_time-1):
        # growth and grazing rates
        u_1next=PO4[t-1]/(kP_1+PO4[t-1])*umax_1
        u_2next=PO4[t-1]/(kP_2+PO4[t-1])*umax_2
       
        u_1.append(u_1next)
        u_2.append(u_2next)
       
        g_1next=P_1[t-1]/(kZ+P_1[t-1]+P_2[t-1])*gmax_1
        g_2next=P_2[t-1]/(kZ+P_1[t-1]+P_2[t-1])*gmax_2
       
        g_1.append(g_1next)
        g_2.append(g_2next)
   
   	    # fluxes
        PP_1next=u_1[t]*P_1[t-1]
        PP_2next=u_2[t]*P_2[t-1]
       
        PP_1.append(PP_1next)
        PP_2.append(PP_2next)
       
        G1_next=g_1[t]*Z[t-1]
        G2_next=g_2[t]*Z[t-1]
        death_Znext=mZ*Z[t-1]**2
        exc_Znext=eZ*Z[t-1]
       
        G1.append(G1_next)
        G2.append(G2_next)
        death_Z.append(death_Znext)
        exc_Z.append(exc_Znext)
       
        PO4_regnext=epsilon*exc_Z[t]
        p_fecnext=(1-epsilon)*exc_Z[t]
       
        PO4_reg.append(PO4_regnext)
        p_fec.append(p_fecnext)
       
       
        #PP_1 puise en premier dans le pool de nutriment
        max_available_PO4 = PO4[t-1]+Psupply[t]*dt+PO4_reg[t]*dt
        if PP_1[t]>max_available_PO4/dt:
            PP_1[t]=max_available_PO4/dt
            
        max_available_PO4=max_available_PO4-PP_1[t]*dt
        if PP_2[t]>max_available_PO4/dt:
            PP_2[t]=max_available_PO4/dt
        
        # export 
        Export_next =p_fec[t]+death_Z[t]
        Export.append(Export_next)
        
        # carbon-based nutrients and biomass
        PO4_next=PO4[t-1]+Psupply[t]*dt+PO4_reg[t]*dt-PP_1[t]*dt-PP_2[t]*dt
        if PO4_next<=0:
            PO4_next=0
            
        P_1_next=P_1[t-1]+PP_1[t]*dt-G1[t]*dt
        if P_1_next<=0:
            P_1_next=0
        
        P_2_next=P_2[t-1]+PP_2[t]*dt-G2[t]*dt
        if P_2_next<=0:
            P_2_next=0
            
        Z_next=Z[t-1]+G1[t]*dt+G2[t]*dt-exc_Z[t]*dt-death_Z[t]*dt
        if Z_next<=0:
            Z_next=0
        
        P_1.append(P_1_next)
        P_2.append(P_2_next)
        Z.append(Z_next)
        PO4.append(PO4_next)
    
    
    # Figures
    plt.figure(1)
    plt.plot(time, P_1, label='P_1', color="chartreuse")
    plt.plot(time, P_2, label='P_2',color="green")
    plt.plot(time, Z, label='Zoo', color="aqua")
    plt.plot(time,PO4, label='PO4',color="magenta")
    plt.xlabel('Time')
    plt.ylabel('mmolC.m^-3')
    plt.legend()
    plt.title('model outputs (concentration over time)')
    
    
    plt.figure(2)
    # Calcul des dérivées de P_1 et P_2 par rapport au temps
    dp1_dt = np.gradient(P_1, dt)
    dp2_dt = np.gradient(P_2, dt)
    
    # Tracé du portrait de phase
    plt.plot(P_2, P_1)
    plt.xlabel('P_2')
    plt.ylabel('P_1')
    
    # Tracé du champ de vecteur sur le portrait de phase
    plt.quiver(P_2, P_1, dp2_dt, dp1_dt)
    

    
    # plt.figure(2)
    # plt.plot(time,Export,label='export',color="black")
    # plt.plot(time,Psupply,label='Psupply',color="red")
    # plt.xlabel('Time (d)')
    # plt.ylabel('Flux (mmolC.m-3.d-1)')
    # plt.legend()
    # plt.title('Budget')
    
    # plt.figure(3)
    # plt.plot(u_1,P_1,label='P_1')
    # plt.plot(u_2,P_2, label='P_2')
    # plt.xlabel('Growth rate')
    # plt.ylabel('Biomass')
    # plt.legend()
    # plt.title('Behrenfeld model outputs')
    
    
    
#Call the function       
ga_model_2P1Z_v8()