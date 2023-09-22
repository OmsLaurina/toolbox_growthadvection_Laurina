#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 15:37:57 2023

@author: loms
"""
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science','no-latex'])
plt.close('all')

# Load the std_variations_abs data from the text files
test_values = ['0.01', '0.1']
std_variations = []
coeff_variations = []

for test in test_values:
    std_variations.append(np.loadtxt(f'std_variations_{test}.txt'))
    coeff_variations.append(np.loadtxt(f'coeff_variations_{test}.txt'))

############ Absolute Standard Deviation 

plt.rc('font', size=10) 
fig, axs = plt.subplots(1, 2, figsize=(10, 5))
plt.subplots_adjust(wspace=0, hspace=0)

param_names2 = [r'$\mu^{max}_1$', r'$\mu^{max}_2$', r'$K^P_1$', r'$K^P_2$', r'$g^{max}_1$', r'$g^{max}_2$', r'$K^{Z}_1$', r'$K^{Z}_2$']
for i, test in enumerate(test_values):
    
    std_variations_abs = std_variations[i]
    cax = axs[i].matshow(std_variations_abs, cmap='plasma')
    axs[i].set_xticks(np.arange(5))
    axs[i].set_yticks(np.arange(len(param_names2)))
    axs[i].set_xticklabels([r'$P_1$', r'$P_2$', r'$Z$', r'$PO4$', r'$R$'])
    axs[i].set_yticklabels(param_names2)
    axs[i].xaxis.tick_bottom()
    axs[i].set_xlabel('State variable')
    axs[i].set_ylabel('Parameter')
    axs[i].set_title(rf'$\mathbf{{P_{{supply}} = {test}}}$', fontweight='bold', fontsize=12)

cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.7])  # [x, y, width, height]
cbar = plt.colorbar(cax, cax=cbar_ax, label='Absolute Standard Deviation')
cbar.set_label('Absolute standard variation', fontsize=10)
cax.set_clim(0, 0.1)
plt.tight_layout()
plt.savefig("../outputs/std.pdf", format='pdf')

############# Coefficient of variation

plt.rc('font', size=10) 
fig, axs = plt.subplots(1, 2, figsize=(10, 5))
plt.subplots_adjust(wspace=0, hspace=0)

param_names2 = [r'$\mu^{max}_1$', r'$\mu^{max}_2$', r'$K^P_1$', r'$K^P_2$', r'$g^{max}_1$', r'$g^{max}_2$', r'$K^Z1$', r'$K^Z2$']
for i, test in enumerate(test_values):
    
    coeff_variations_abs = coeff_variations[i]
    cax = axs[i].matshow(coeff_variations_abs, cmap='plasma')
    axs[i].set_xticks(np.arange(5))
    axs[i].set_yticks(np.arange(len(param_names2)))
    axs[i].set_xticklabels([r'$P_1$', r'$P_2$', r'$Z$', r'$PO4$', r'$R$'])
    axs[i].set_yticklabels(param_names2)
    axs[i].xaxis.tick_bottom()
    axs[i].set_xlabel('State variable')
    axs[i].set_ylabel('Parameter')
    axs[i].set_title(rf'$\mathbf{{P_{{supply}} = {test}}}$', fontweight='bold', fontsize=12)

cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.7])  # [x, y, width, height]
cbar = plt.colorbar(cax, cax=cbar_ax, label='Absolute Standard Deviation')
cbar.set_label('Variation coefficient', fontsize=10)
cax.set_clim(0, 10)
plt.tight_layout()
plt.savefig("../outputs/coefvar.pdf", format='pdf')