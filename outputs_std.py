#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 15:32:50 2023

@author: loms
"""
import numpy as np
import numpy as npy
import matplotlib.pyplot as plt
from growth_model_2P1Z_v10 import growth_model_2P1Z_v10
import scienceplots
plt.style.use(['science','no-latex'])
plt.close('all')

dt = 0.2
end_time = 2000
time = npy.arange(0, end_time, dt)
p = 0.01
test=p
Psupply = [p] * len(time)

param_names = ['umax1', 'umax2', 'kP1', 'kP2', 'gmax1', 'gmax2', 'kZ1', 'kZ2']

# Default values
default_params = {
    'umax1': 1.9872,
    'umax2': 2.7648,
    'kP1': 1,
    'kP2': 3,
    'gmax1': 3.89,
    'gmax2': 0.43,
    'kZ1': 5,
    'kZ2': 20,
}


tx = 0.5  # Pourcentage of variation around the default value
n = 5 # Number of value tested

param_ranges = [
    np.linspace(default_params['umax1'] * (1 - tx), default_params['umax1'] * (1 + tx), n),  # umax1
    np.linspace(default_params['umax2'] * (1 - tx), default_params['umax2'] * (1 + tx), n),  # umax2
    np.linspace(default_params['kP1'] * (1 - tx), default_params['kP1'] * (1 + tx), n),      # kP1
    np.linspace(default_params['kP2'] * (1 - tx), default_params['kP2'] * (1 + tx), n),      # kP2
    np.linspace(default_params['gmax1'] * (1 - tx), default_params['gmax1'] * (1 + tx), n),  # gmax1
    np.linspace(default_params['gmax2'] * (1 - tx), default_params['gmax2'] * (1 + tx), n),  # gmax2
    np.linspace(default_params['kZ1'] * (1 - tx), default_params['kZ1'] * (1 + tx), n),      # kZ1
    np.linspace(default_params['kZ2'] * (1 - tx), default_params['kZ2'] * (1 + tx), n),      # kZ2
]

# Get the default solutions
P1_default, P2_default, Z_default, PO4_default, _ = growth_model_2P1Z_v10(Psupply, time)
P1_default = P1_default[-1]
P2_default = P2_default[-1]
Z_default = Z_default[-1]
PO4_default = PO4_default[-1]
R_default = P1_default / (P1_default + P2_default)

# Initialize empty matrices to store the std and the variation coefficient
variations = np.zeros((5,len(param_names),n)) # 5 state variables

fig, axs = plt.subplots(5, 1, figsize=(10, 15))

custom_colors =['#FFB6C1', '#FFD700', '#98FB98', '#87CEEB', '#FFA07A', '#9370DB', '#90EE90', '#F0E68C']

for p_idx, p_name in enumerate(param_names):
    
    # Loop along the parameter's values
    for p_val_idx, p_value in enumerate(param_ranges[p_idx]):
        
        # Copy parameters and modify the value at each iteration
        modified_params = default_params.copy()
        modified_params[p_name] = p_value

        # Get the solutions with the new parameter's value
        P1, P2, Z, PO4, _ = growth_model_2P1Z_v10(Psupply, time, **modified_params)
        
        # Plot P1
        axs[0].plot(time, P1, label=f'{p_name}={p_value}',color=custom_colors[p_idx])
        axs[0].axhline(y=P1_default, color='red', linestyle='--', label=f'Default {p_name}')
        axs[0].set_title('P1')
        
        # Plot P2
        axs[1].plot(time, P2, label=f'{p_name}={p_value}',color=custom_colors[p_idx])
        axs[1].axhline(y=P2_default, color='red', linestyle='--', label=f'Default {p_name}')
        axs[1].set_title('P2')
        
        # Plot R
        R = np.array(P1) / (np.array(P1) + np.array(P2))
        axs[2].plot(time, R, label=f'{p_name}={p_value}',color=custom_colors[p_idx])
        axs[2].axhline(y=R_default, color='red', linestyle='--', label=f'Default {p_name}')
        axs[2].set_title('R')
        
        # Plot Z
        axs[3].plot(time, Z, label=f'{p_name}={p_value}',color=custom_colors[p_idx])
        axs[3].axhline(y=Z_default, color='red', linestyle='--', label=f'Default {p_name}')
        axs[3].set_title('Z')
        
        # Plot PO4
        axs[4].plot(time, PO4, label=f'{p_name}={p_value}',color=custom_colors[p_idx])
        axs[4].axhline(y=PO4_default, color='red', linestyle='--', label=f'Default {p_name}')
        axs[4].set_title('PO4')


        # Calculate the absolute standard deviation and the variation coefficient 
        
        P1 = P1[-1]
        P2 = P2[-1]
        Z = Z[-1]
        PO4 = PO4[-1]
        R = P1/(P1 + P2)

        variations[:,p_idx, p_val_idx] = np.array([
            (P1 - P1_default),
            (P2 - P2_default),
            (Z - Z_default),
            (PO4 - PO4_default),
            (R - R_default)
        ])
        
        std_variations = np.abs(np.std(variations, axis=2))
        std_variations = std_variations.T
        coeff_variations = np.abs(np.std(variations, axis=2))/np.abs(np.mean(variations, axis=2))
        coeff_variations = coeff_variations.T

plt.subplots_adjust(hspace=0.5)
plt.show()
plt.savefig(f'../outputs/sensitivity_test_{p}.pdf', format='pdf')

# Save data to a text file
np.savetxt(f'../outputs/std_variations_{test}.txt', std_variations, fmt='%.8f')
np.savetxt(f'../outputs/coeff_variations_{test}.txt', coeff_variations, fmt='%.8f')
