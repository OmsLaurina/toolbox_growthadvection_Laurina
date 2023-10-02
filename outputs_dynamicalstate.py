#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 16:51:23 2023

@author: loms
"""

from growth_model_2P1Z_v10 import growth_model_2P1Z_v10 
from fluxes import periodicflux, pulsedfluxGrover1990, pulsedflux_randomAmplitude,pulsedflux_stepfunction
import matplotlib.pyplot as plt
import numpy as npy
import numpy as np

""" Set up """

# Timestep (same than timestep of the growth function : by default is 0.2)
dt = 0.1
# Choose the time at which the simulation end (days) (steady-state)
end_time = 2000
# create a time vector from 0 to end_time, seperated by a timestep
time = npy.arange(1,end_time,dt)


""" FLUXES """
# Externe input of nutrients 
Psupply_ini = 0.05 #1 #0.04419 = valeur ou R = 50/50
#Psupply variable
moy_flux = Psupply_ini #1 #0.01
Psupply_cst =[Psupply_ini] * len(time)
#Amplitude
b = 0.02
#Period
T = 5/dt # Create a figure with subplots
nbpulse = 1
   
#Time of the simulated periodic fluxes
end_time_flux = 1999
time_flux = npy.arange(1,end_time_flux,dt)
nb_time=len(time_flux)

fluxes = 'pulsedflux_stepfunction'

#[Psupply,arg] = periodicflux(moy_flux, time_flux, b=b, T=T)
#[Psupply,arg] = pulsedfluxGrover1990(moy_flux, time_flux)
#[Psupply,arg] = pulsedflux_randomAmplitude(moy_flux, time_flux, T=T)

if nbpulse == 3 : 
    [Psupply,arg] = pulsedflux_stepfunction(time_flux, A=b,C=moy_flux, pulsations=[{'t1': 5, 't2': 10},{'t1': 20, 't2': 25},{'t1': 35, 't2': 40}])
if nbpulse == 2 : 
    [Psupply,arg] = pulsedflux_stepfunction(time_flux, A=b,C=moy_flux, pulsations=[{'t1': 5, 't2': 10},{'t1': 20, 't2': 25}])
if nbpulse == 1 : 
    [Psupply,arg] = pulsedflux_stepfunction(time_flux, A=b,C=moy_flux, pulsations=[{'t1': 5, 't2': 10}])
    
    
""" Solutions in steady state """
# Get the solution ofthe model at the steady state equilibrium 
[P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply_cst,time,dt)
P1_new = P1[len(P1)-1]
P2_new = P2[len(P2)-1]
Z_new = Z[len(Z)-1]
PO4_new = PO4[len(PO4)-1]
ratio = P1[len(P1)-1]/(P1[len(P1)-1]+P2[len(P2)-1])

""" Solutions in dynamical state """
# Run the model from the equilibrium solutions with the dynamical fluxes of nutrient
[P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply,time_flux,dt,P1_ini=P1_new,P2_ini=P2_new,Z_ini=Z_new,PO4_ini=PO4_new)

Psupply_list = Psupply.tolist()

data_filename = f'../outputs/data_pulse{nbpulse}_b{b}_d{end_time_flux}.txt'
data = np.column_stack((time_flux, P1, P2, Z, PO4,Psupply_list))
np.savetxt(data_filename, data, header="time_flux P1 P2 Z PO4 Psupply")