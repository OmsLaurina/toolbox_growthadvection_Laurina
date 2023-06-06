"""
Test on variable input of nutrient
"""

import matplotlib.pyplot as plt
import numpy as npy
import numpy as np
import scienceplots
from matplotlib.cm import ScalarMappable
import sys
sys.path.append('../')

plt.style.use(['science','no-latex'])
plt.close('all')

from growth_model_2P1Z_v10 import growth_model_2P1Z_v10 

""" Set up """

# Timestep (same than timestep of the growth function : by default is 0.2)
dt = 0.2
# Choose the time at which the simulation end (days) (steady-state)
end_time = 100
# create a time vector from 0 to end_time, seperated by a timestep
time = npy.arange(0,end_time,dt)

# Externe input of nutrients 
Psupply_moy = 0.01 #1 #0.01

#Psupply constant
Psupply_cst =[Psupply_moy] * len(time)

#Psupply periodic
Psupply = []
b = 1 #0.01
T = 30
nb_time=len(time)
for i in range(0,nb_time):
      Psupplynext = Psupply_moy+b*np.sin(2*np.pi*i/T)
      Psupply.append(Psupplynext)

 
P1_new=np.zeros((len(Psupply)))
P2_new=np.zeros((len(Psupply)))
Z_new=np.zeros((len(Psupply)))
PO4_new=np.zeros((len(Psupply))) 
ratio=np.zeros((len(Psupply)))

[P1,P2,Z,PO4,Export,u1,u2,g1,g2,arg]=growth_model_2P1Z_v10(Psupply_cst,time)
P1_new = P1[len(P1)-1]
P2_new = P2[len(P2)-1]
Z_new = Z[len(Z)-1]
PO4_new = PO4[len(PO4)-1]
ratio = P1[len(P1)-1]/(P1[len(P1)-1]+P2[len(P2)-1])

[P1,P2,Z,PO4,Export,u1,u2,g1,g2,arg]=growth_model_2P1Z_v10(Psupply,time,P1_ini=P1_new,P2_ini=P2_new,Z_ini=Z_new,PO4_ini=PO4_new)

    
plt.figure(1)

plt.plot(time,Psupply, label='Psupp', color="blue")
plt.plot(time, P1, label='P1', color="chartreuse")
plt.plot(time, P2, label='P2', color="green")
plt.plot(time, Z, label='Z', color="aqua")
plt.plot(time,PO4, label='PO4', color="magenta")
plt.xlabel('Time [day]')
plt.ylabel('Biomasses [$mmolC.m^{-3}$]')
plt.legend(frameon=True) 