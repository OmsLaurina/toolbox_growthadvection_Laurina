#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 11:40:23 2023

@author: loms
"""
import numpy as npy
import numpy as np
from growth_model_2P1Z_v10 import growth_model_2P1Z_v10

""" Set up """

grazing = 'diffgrazing'

# Timestep (same than timestep of the growth function : by default is 0.2)
dt = 0.1
end_time = 2000
# create a time vector from 0 to end_time, seperated by a timestep
time = npy.arange(0,end_time,dt)
nb_time = len(time)

# Externe input of nutrients 
#!! WARNING DEFINE PSUPPLY = 1 IF YOU WANT TO TEST THIS PARAMETER !!
Psupply_moy = 1 #1 #0.01
Psupply =[Psupply_moy] * len(time)
# !!!! use Psupply_arr instead of Psupply when you call the function to test a range of Psupply values (Psupply_moy must be equal to 1)!!!
Psupply_arr = np.array(Psupply)

""" MODIFY THIS TO CHOOSE THE PARAMETERS AND VALUE RANGE YOU WANT TO TEST """

number_param = 1 #1,2
n = 100 #number of value tested

#name of the parameters in the model you want to study 
param1 = 'Psupp'
param2 = 'gmax'

# name for plot
name_param1=r'$P_{supply}$'
name_param2=r'$g_{max}$'

#Tested range of values
min_param1 = 0.01
max_param1 = 0.1
l_param1 = np.linspace(min_param1,max_param1,n)

min_param2 = 0.1
max_param2 = 5
l_param2 = np.linspace(min_param2,max_param2,n)

if number_param == 2:
    
    filename=f"../outputs/data_sensitivitytest_{param1}_{param2}.txt"
    
    """ loop on i for param1 and j for param2 """ 
    
    P_1=np.zeros((len(l_param1),len(l_param2)))
    P_2=np.zeros((len(l_param1),len(l_param2)))
    ratio=np.zeros((len(l_param1),len(l_param2)))
    
    # !!!! WARNING : Modify the kwargs name to select the good parameters !!!! #

    i = 0
    for param1 in l_param1:
        j=0
        for param2 in l_param2:
            [P1,P2,Z,PO4,Export,arg]=growth_model_2P1Z_v10(Psupply_arr,time, dt, gmax1=param2,gmax2=param2,kZ1=param1,kZ2=param1) #A modifier
            P_1[i,j] = P1[len(P1)-1]
            P_2[i,j] = P2[len(P2)-1]
            ratio[i,j] = P1[len(P1)-1]/(P1[len(P1)-1]+P2[len(P2)-1])
            j=j+1
        i=i+1
        print(i)
    
else:
    filename=f"../outputs/data_senstitivitytest_{param1}_{grazing}.txt"
    
    'If you want to study only param1'
    
    P_1=np.zeros((len(l_param1)))
    P_2=np.zeros((len(l_param1)))
    ratio=np.zeros((len(l_param1)))
    P_O4=np.zeros((len(l_param1)))
    prop_P1=np.zeros((len(l_param1)))
    prop_P2=np.zeros((len(l_param1)))
    
    # !!!! WARNING : Modify the kwargs name to select the good parameters !!!! #     
    
    i=0
    for param1 in l_param1:
        if grazing == 'diffgrazing':
            [P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply_arr*param1,time, dt)
        if grazing == 'equalgrazing':
            [P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply_arr*param1,time, dt, gmax2=3.89,kZ2=5)#A modifier
        if grazing == 'nograzing':
            [P1,P2,Z,PO4,arg]=growth_model_2P1Z_v10(Psupply_arr*param1,time, dt, gmax1=0,gmax2=0)#A modifier
            
        P_1[i] = P1[len(P1)-1]
        P_2[i] = P2[len(P2)-1]
        ratio[i] = P1[len(P1)-1]/(P1[len(P1)-1]+P2[len(P2)-1])
        P_O4[i] = PO4[len(PO4)-1]
        prop_P1[i] = (P1[len(P1)-1]/(P1[len(P1)-1]+P2[len(P2)-1]))*100
        prop_P2[i] =  (P2[len(P2)-1]/(P1[len(P1)-1]+P2[len(P2)-1]))*100
        i=i+1
        print(i)

# Pour calculer les proportions
for index, param1_value in enumerate(l_param1):
    if param1_value < 0.0509:
        prop_p1_value = prop_P1[index]
        prop_p2_value = prop_P2[index]
        
        with open("prop_p1_param1_inf_0.05.txt", "a") as f_prop_p1_inf_005:
            f_prop_p1_inf_005.write(f"{prop_p1_value}\t{param1_value}\n")  # Écrit la valeur de param1 en deuxième colonne
        with open("prop_p2_param1_inf_0.05.txt", "a") as f_prop_p2_inf_005:
            f_prop_p2_inf_005.write(f"{prop_p2_value}\t{param1_value}\n")  # Écrit la valeur de param1 en deuxième colonne
        print(index)
        print(param1_value)
    else:
        prop_p1_value = prop_P1[index]
        prop_p2_value = prop_P2[index]
        
        with open("prop_p1_param1_sup_0.05.txt", "a") as f_prop_p1_sup_005:
            f_prop_p1_sup_005.write(f"{prop_p1_value}\t{param1_value}\n")  # Écrit la valeur de param1 en deuxième colonne
        with open("prop_p2_param1_sup_0.05.txt", "a") as f_prop_p2_sup_005:
            f_prop_p2_sup_005.write(f"{prop_p2_value}\t{param1_value}\n")  # Écrit la valeur de param1 en deuxième colonne
        print(index)
        print(param1_value)
 
filename_P1 = filename.replace('.txt', '_P1.txt')
filename_P2 = filename.replace('.txt', '_P2.txt')
filename_ratio = filename.replace('.txt', '_ratio.txt')
filename_PO4 = filename.replace('.txt', '_PO4.txt')

np.savetxt(filename_P1, P_1)
np.savetxt(filename_P2, P_2)
np.savetxt(filename_ratio, ratio)
np.savetxt(filename_PO4, P_O4)