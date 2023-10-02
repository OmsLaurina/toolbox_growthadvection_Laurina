#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 11:18:36 2023

@author: loms
"""

import numpy as np
from growth_model_2P1Z_v10 import growth_model_2P1Z_v10

# Set up
dt = 0.1
end_time = 2000
time = np.arange(0, end_time, dt)
Psupply_moy = 0.01
Psupply = [Psupply_moy] * len(time)

test = Psupply_moy

# Call the function
[P1, P2, Z, PO4, arg] = growth_model_2P1Z_v10(Psupply, time)

# # Save data to a text file
# data_filename = f'../outputs/data_{test}.txt'
# data = np.column_stack((time, P1, P2, Z, PO4))
# np.savetxt(data_filename, data, header="Time P1 P2 Z PO4")

# ##### FOR SCENARIOS #####

# grazing = 'nograzing'

# [P1, P2, Z, PO4, arg] = growth_model_2P1Z_v10(Psupply, time, gmax1=0, gmax2=0)
# data_filename = f'../outputs/data_{grazing}_{test}.txt'
# data = np.column_stack((time, P1, P2, Z, PO4))
# # np.savetxt(data_filename, data, header="Time P1 P2 Z PO4")

# grazing = 'equalgrazing'

# [P1, P2, Z, PO4, arg] = growth_model_2P1Z_v10(Psupply, time, kZ2=5, gmax2=3.89)
# data_filename = f'../outputs/data_{grazing}_{test}.txt'
# data = np.column_stack((time, P1, P2, Z, PO4))
# # np.savetxt(data_filename, data, header="Time P1 P2 Z PO4")

grazing = 'diffgrazing'

[P1, P2, Z, PO4, arg] = growth_model_2P1Z_v10(Psupply, time, dt)
prop_P1 = (P1[-1]/(P1[-1]+P2[-1]))*100
prop_P2 =  (P2[-1]/(P1[-1]+P2[-1]))*100

data_filename = f'../outputs/data_{grazing}_{test}.txt'
data = np.column_stack((time, P1, P2, Z, PO4))
np.savetxt(data_filename, data, header="Time P1 P2 Z PO4")

