#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 17:24:51 2023

@author: loms
"""
import numpy as npy
import numpy as np
from growth_model_2P1Z_v10 import growth_model_2P1Z_v10 
from fluxes import periodicflux, pulsedfluxGrover1990, pulsedflux_randomAmplitude,pulsedflux_stepfunction


"""  ###### Set up ###### """

# Timestep (same than timestep of the growth function : by default is 0.2)
dt = 0.1

#Name test
test="pulsedflux_stepfunction" #pulsed #periodic

# Externe input of nutrients 
Psupply_ini = 0.05#for the steady state model, 0.04419
moy_flux = 0.05#mean of the sinusoïdale function

#Amplitude
b = 0.08 #0.01
#Period
T = 5/dt

#Time of the simulated periodic fluxes
end_time_flux = 90
time_flux = npy.arange(0,end_time_flux,dt)
nb_time=len(time_flux)
nbpulse = 1

"""  ###### END SET UP ###### """
##################################################################################################################
"""  ###### Get the solution of the model at the steady state equilibrium ###### """

# Choose the time at which the simulation end (days) (steady-state)
end_time = 2000
# create a time vector from 0 to end_time, seperated by a timestep
time = npy.arange(0,end_time,dt)

#Psupply constant
Psupply_cst =[Psupply_ini] * len(time)
[P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply_cst,time, dt)
P1_eq = P1[len(P1)-1]
P2_eq = P2[len(P2)-1]
Z_eq = Z[len(Z)-1]
PO4_eq = PO4[len(PO4)-1]

"""  ###### END Get the solution of the model at the steady state equilibrium ###### """
##################################################################################################################
"""  ###### Run the model from the equilibrium solutions with the dynamical fluxes of nutrient ###### """

number_param = 1

n = 10 #number of value tested
n2 = n #to reduce the number of line in the line plots 

#name of the parameters in the model you want to study 
#param1 = r'$\overline{P_{supply}}$' 
param1_bis = 'Psupply_moy'
param1='b'
param2 = 'T'

# name for plot
name_param1=param1
name_param2=param2

#Tested range of values
min_param1 = 0.01
max_param1 = 0.1
l_param1 = np.linspace(min_param1,max_param1,n)
    
min_param2 = 1
max_param2 = 30
l_param2 = np.linspace(min_param2,max_param2,n)


# !!!! WARNING : Modify the kwargs name to select the good parameters !!!! #

if number_param == 2:
    
    filename=f"../outputs/data_senstitivitytest_{param1}_{param2}pulse{nbpulse}_b{b}_d{end_time_flux}.txt"
    
    """ loop on i for param1 and j for param2 """ 
    
    P1_new = np.zeros((len(l_param1), len(l_param2), nb_time))
    P2_new = np.zeros((len(l_param1), len(l_param2), nb_time))
    ratio = np.zeros((len(l_param1), len(l_param2), nb_time))
    Psupply_new = np.zeros(nb_time)
    
    i = 0
    for param1 in l_param1:
        j = 0
        for param2 in l_param2:
            
            if test=='periodicflux':
                [Psupply,arg] = periodicflux(moy_flux, time_flux, b=param1, T=param2)
                [P1,P2,Z,PO4,Export,u1,u2,g1,g2,arg]=growth_model_2P1Z_v10(Psupply,time_flux,P1_ini=P1_eq,P2_ini=P2_eq,Z_ini=Z_eq,PO4_ini=PO4_eq)  
            elif test=='pulsedfluxGrover1990':
                [Psupply,arg] = pulsedfluxGrover1990(moy_flux, time_flux)
                [P1,P2,Z,PO4,Export,u1,u2,g1,g2,arg]=growth_model_2P1Z_v10(Psupply,time_flux,P1_ini=P1_eq,P2_ini=P2_eq,Z_ini=Z_eq,PO4_ini=PO4_eq) 
            elif test=='pulsedflux_randomAmplitude':
                [Psupply,arg] = pulsedflux_randomAmplitude(moy_flux, time_flux, T=param2)
                [P1,P2,Z,PO4,Export,u1,u2,g1,g2,arg]=growth_model_2P1Z_v10(Psupply,time_flux,P1_ini=P1_eq,P2_ini=P2_eq,Z_ini=Z_eq,PO4_ini=PO4_eq)
            elif test=='pulsedflux_stepfunction':
                 [Psupply,arg] = pulsedflux_stepfunction(time_flux, A=param1,C=moy_flux, pulsations=[{'t1': 5, 't2': 10}])
                 [P1,P2,Z,PO4,Export,u1,u2,g1,g2,arg]=growth_model_2P1Z_v10(Psupply,time_flux, dt, P1_ini=P1_eq,P2_ini=P2_eq,Z_ini=Z_eq,PO4_ini=PO4_eq)
                
            P1_new[i, j, :] = np.array(P1)
            P2_new[i, j, :] = np.array(P2)
            ratio[i, j, :] =np.array(P1) / (np.array(P1) + np.array(P2))
        
            j += 1
        i += 1
        print(i)
    sum_ratios = np.sum(ratio, axis=2)
    mean_ratios = sum_ratios/nb_time
    
else:
    filename=f"../outputs/data_senstitivitytest_{param1}_pulse{nbpulse}_b{b}_d{end_time_flux}.txt"
    
    'If you want to study only param1'
    pulsedflux_stepfunction
    P1_new = np.zeros((len(l_param1), nb_time))
    P2_new = np.zeros((len(l_param1), nb_time))
    PO4_new = np.zeros((len(l_param1), nb_time))
    Z_new = np.zeros((len(l_param1), nb_time))
    ratio = np.zeros((len(l_param1),nb_time))
    Psupply_new =  np.zeros((len(l_param1), nb_time))
    b_list=np.zeros(len(l_param1))
     
    i=0
    for param1 in l_param1:
        
       if test=='periodicflux':
           [Psupply,arg] = periodicflux(moy_flux*param1, time_flux, b=b, T=T)
           [P1,P2,Z,PO4,Export,u1,u2,g1,g2,arg]=growth_model_2P1Z_v10(Psupply,time_flux,P1_ini=P1_eq,P2_ini=P2_eq,Z_ini=Z_eq,PO4_ini=PO4_eq) 
       elif test=='pulsedfluxGrover1990':
           [Psupply,arg] = pulsedfluxGrover1990(moy_flux*param1, time_flux)
           [P1,P2,Z,PO4,Export,u1,u2,g1,g2,arg]=growth_model_2P1Z_v10(Psupply,time_flux,P1_ini=P1_eq,P2_ini=P2_eq,Z_ini=Z_eq,PO4_ini=PO4_eq)   
       elif test=='pulsedflux_randomAmplitude':
            [Psupply,arg] =  pulsedflux_randomAmplitude(moy_flux*param1, time_flux, T=T)
            [P1,P2,Z,PO4,Export,u1,u2,g1,g2,arg]=growth_model_2P1Z_v10(Psupply,time_flux,P1_ini=P1_eq,P2_ini=P2_eq,Z_ini=Z_eq,PO4_ini=PO4_eq)
       elif test=='pulsedflux_stepfunction':
           if nbpulse==1:
            [Psupply,arg] = pulsedflux_stepfunction(time_flux, A=param1,C=moy_flux, pulsations=[{'t1': 5, 't2': 10}])
            [P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply,time_flux, dt, P1_ini=P1_eq,P2_ini=P2_eq,Z_ini=Z_eq,PO4_ini=PO4_eq)
           if nbpulse ==2:
               [Psupply,arg] = pulsedflux_stepfunction(time_flux, A=param1,C=moy_flux, pulsations=[{'t1': 5, 't2': 10},{'t1': 20, 't2': 30}])
               [P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply,time_flux, dt, P1_ini=P1_eq,P2_ini=P2_eq,Z_ini=Z_eq,PO4_ini=PO4_eq)
           if nbpulse ==3:
                [Psupply,arg] = pulsedflux_stepfunction(time_flux, A=param1,C=moy_flux, pulsations=[{'t1': 5, 't2': 10},{'t1': 20, 't2': 30},{'t1': 40, 't2': 50}])
                [P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply,time_flux, dt, P1_ini=P1_eq,P2_ini=P2_eq,Z_ini=Z_eq,PO4_ini=PO4_eq)
      
       P1_new[i, :] = np.array(P1)
       P2_new[i, :] = np.array(P2)
       PO4_new[i, :] = np.array(PO4)
       Z_new[i, :] = np.array(Z)
       ratio[i,:] = np.array(P1)/(np.array(P1)+ np.array(P2))
       Psupply_new[i, :] = np.array(Psupply)
       i=i+1
       print(i)
       
filename_P1 = filename.replace('.txt', '_P1.txt')
filename_P2 = filename.replace('.txt', '_P2.txt')
filename_ratio = filename.replace('.txt', '_ratio.txt')
filename_PO4 = filename.replace('.txt', '_PO4.txt')
filename_Z = filename.replace('.txt', '_Z.txt')
filename_Psupply = filename.replace('.txt', '_Psupply.txt')
        
np.savetxt(filename_P1, P1_new)
np.savetxt(filename_P2, P2_new)
np.savetxt(filename_ratio, ratio)
np.savetxt(filename_PO4, PO4_new)
np.savetxt(filename_Z, Z_new)
np.savetxt(filename_Psupply, Psupply_new)