#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Implementation of the Behrenfeld Model 
                2P1Z
with growth rate = monod equation
                
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
u_big=[]			
u_small=[]	
#Primary production
PP_small=[]
PP_big=[]
# biomass of Phytoplankton and Zooplankton (mmolC/m³)
P_small=[] 			
P_big=[] 						
Z = []
#Phosphate concentration
PO4 = []
#Excretion and regeneration
exc_Z=[]
PO4_reg=[]
# initial conditions at time=0	
u0_small = 0
u0_big = 0
P0_small = 0.6
P0_big = 0.1
Z0 = 0.3
PO40 = 0.5

"""" definition of model parameters"""
##Berhenfeld parameters
# Herbivore grazing rate (d^-1)
c1_small = 0.1
c1_big= 0.2
# Ingestion efficiency
c2 = 0.2
# Predatory loss rate (d^-1)
c3  = 0.005

##Monod parameters
kP_small = 13.3
kP_big = 15.3
umax_small = 1.9872
umax_big = 2.4226 

##Growth model parameters
#zoo excretion fraction
eZ = 0.1
#fraction of Zbig excretion that is available as regenerated PO4
epsilon = 0.75

# lag-time between division and loss rates
j = 1

#Apport de nutriment externe
Psupply = np.linspace(0.001,0.1,len(time))


def Behrenfeld_model_1P1Z_v1():
    
    ## -------------- Parameters
    global dt, end_time, time, c1_small, c1_big, c2, c3, u_small, u_big, P_small, P_big, Z, u0_small, u0_big, P0_small, P0_big, Z0, j, Psupply
    
    # Initial conditions at time=0	
    u_small.append(u0_small)
    u_big.append(u0_big)		
    P_small.append(P0_small)
    P_big.append(P0_big)
    Z.append(Z0)
    PO4.append(PO40)

    ## -------------- loop on time

    nb_time=len(time)
    
    for t in range(1,nb_time):
        
        #Primary production (=uptake nutrients)
        PP_small_next = u_small[t-1]*P_small[t-1]
        PP_big_next = u_big[t-1]*P_big[t-1]
        
        PP_small.append(PP_small_next)
        PP_big.append(PP_big_next)
        
        #Phytoplankton biomass
        P_small_next = P_small[t-1] + PP_small[t-1]*dt - (c1_small*P_small[t-1]*Z[t-1])*dt
        P_big_next = P_big[t-1] + PP_big[t-1]*dt - (c1_big*P_big[t-1]*Z[t-1])*dt
        
        #Excretion and regenration
        exc_Z_next =  eZ*Z[t-1]
        PO4_reg_next = epsilon*Z[t-1]
        
        exc_Z.append(exc_Z_next)
        PO4_reg.append(PO4_reg_next)

        #Zooplankton biomass
        Z_next = Z[t-1] + (c1_small*c2*P_small[t-1]*Z[t-1])*dt + (c1_big*c2*P_big[t-1]*Z[t-1])*dt - (c3*Z[t-1]**2)*dt-exc_Z[t-1]*dt
        
        
        #Phosphate concentration
        PO4_next=PO4[t-1]+Psupply[t-1]*dt+PO4_reg[t-1]*dt-PP_small[t-1]*dt-PP_big[t-1]*dt;
        
        P_small.append(P_small_next)
        P_big.append(P_big_next)
        Z.append(Z_next)
        PO4.append(PO4_next)
        
        # #Growth rate calculated with j (ie: lag-time between division and loss rates)
        # if t % j == 0:
        #     u_small_next = PO4[t-j]/(PO4[t-j]+kP_small)*umax_small
        #     u_big_next = PO4[t-j]/(PO4[t-j]+kP_big)*umax_big
        
        u_small_next = PO4[t-1]/(PO4[t-1]+kP_small)*umax_small
        u_big_next = PO4[t-1]/(PO4[t-1]+kP_big)*umax_big
        
        u_big.append(u_big_next)
        u_small.append(u_small_next)
        
        
    # Figures
    plt.figure(1)
    plt.plot(time, P_small, label='P_small')
    plt.plot(time, P_big, label='P_big')
    plt.plot(time, Z, label='Zoo')
    plt.xlabel('Time')
    plt.legend()
    plt.title('Behrenfeld model outputs')
    
    plt.figure(2)
    plt.plot(time,u_small,label='u_small')
    plt.plot(time,u_big, label='u_big')
    plt.xlabel('Time')
    plt.legend()
    plt.title('Behrenfeld model outputs')
    
    plt.figure(3)
    plt.plot(u_small,P_small,label='P_small')
    plt.plot(u_big,P_big, label='P_big')
    plt.xlabel('Growth rate')
    plt.ylabel('Biomass')
    plt.legend()
    plt.title('Behrenfeld model outputs')
    
#Call the function       
Behrenfeld_model_1P1Z_v1() 
