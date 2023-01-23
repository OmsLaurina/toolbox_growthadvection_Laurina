#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Implementation of the Behrenfeld Model 
                2P1Z
                
Created on Tue Dec 20 12:39:07 2022

@author: loms
"""
import matplotlib.pyplot as plt
import numpy as npy


# timestep
dt = 0.2
# the time at which the simulation ends (days)
end_time = 100
# creates a time vector from 0 to end_time, seperated by a timestep
time = npy.arange(0,end_time,dt)


"""" definition of model variables"""
##Creating empty vectors
# growth rates (d^-1)
u_big=[]			
u_small=[]	
# biomass of Phytoplankton and Zooplankton (mmolC/m³)
P_small=[] 			
P_big=[] 						
Z = []
# initial conditions at time=0	
u0_small = 0
u0_big = 0
P0_small = 0.6
P0_big = 0.1
Z0 = 0.3

"""" definition of model parameters"""
# Herbivore grazing rate (d^-1)
c1_small = 0.1
c1_big= 0.2
# Ingestion efficiency
c2 = 0.2
# Predatory loss rate (d^-1)
c3  = 0.005

# lag-time between division and loss rates
j = 1


def Behrenfeld_model_1P1Z_v1():
    
    ## -------------- Parameters
    global dt, end_time, time, c1_small, c1_big, c2, c3, u_small, u_big, P_small, P_big, Z, u0_small, u0_big, P0_small, P0_big, Z0, j
    
    # Initial conditions at time=0
    u_small.append(u0_small)
    u_big.append(u0_big)		
    P_small.append(P0_small)
    P_big.append(P0_big)
    Z.append(Z0)

    ## -------------- loop on time

    nb_time=len(time)
    
    for t in range(1,nb_time):
        
        #Phytoplankton biomass
        P_small_next = P_small[t-1] + (u_small[t-1]*P_small[t-1])*dt - (c1_small*P_small[t-1]*Z[t-1])*dt
        P_big_next = P_big[t-1] +(u_big[t-1]*P_big[t-1])*dt - (c1_big*P_big[t-1]*Z[t-1])*dt
        
        #Zooplankton biomass
        Z_next = Z[t-1] + (c1_small*c2*P_small[t-1]*Z[t-1])*dt + (c1_big*c2*P_big[t-1]*Z[t-1])*dt - (c3*Z[t-1]**2)*dt
        
        P_small.append(P_small_next)
        P_big.append(P_big_next)
        Z.append(Z_next)
        
        #Growth rate calculated with j (ie: lag-time between division and loss rates)
        if t % j == 0:
            u_small_next = u_small[t-j] + (1/P_small[t-1]) * (P_small[t-1]-P_small[t-j])/dt
            u_big_next = u_big[t-j] + (1/P_big[t-1]) * (P_big[t-1]-P_big[t-j])/dt
        
        u_big.append(u_big_next)
        u_small.append(u_small_next)
        
        
    # Figures
    #time evolution of biomasses
    plt.figure(1)
    plt.plot(time, P_small, label='P_small')
    plt.plot(time, P_big, label='P_big')
    plt.plot(time, Z, label='Zoo')
    plt.xlabel('Time')
    plt.legend()
    plt.title('Behrenfeld model outputs')
    
    #time volution of growth rates
    plt.figure(2)
    plt.plot(time,u_small,label='u_small')
    plt.plot(time,u_big, label='u_big')
    plt.xlabel('Time')
    plt.legend('growth rates')
    plt.title('Behrenfeld model outputs')
    
    #biomasses versus growth rates
    plt.figure(3)
    plt.scatter(u_small,P_small,label='P_small')
    plt.scatter(u_big,P_big, label='P_big')
    plt.xlabel('Growth rate')
    plt.ylabel('Biomass')
    plt.legend()
    plt.title('Behrenfeld model outputs')
    
    
#Call the function       
Behrenfeld_model_1P1Z_v1() 
        
