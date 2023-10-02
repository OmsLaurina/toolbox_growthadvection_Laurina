#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 17:35:01 2023

@author: loms
"""

import matplotlib.pyplot as plt
import numpy as npy
import numpy as np
import scienceplots
from matplotlib.cm import ScalarMappable
from matplotlib.gridspec import GridSpec
import sys
sys.path.append('../')

plt.style.use(['science','no-latex'])
plt.close('all')



####### 30 days simulations ########

# Create a figure with subplots
fig, axs = plt.subplots(2, 3, figsize=(10, 5))

dt = 0.1
end_time_flux = 30
time_flux = np.arange(0, end_time_flux, dt)
nb_time = len(time_flux)
min_param = 0.01
max_param = 0.1
l_param = np.linspace(min_param, max_param, 10)

# Load the generated data for the first subplot
data1 = np.loadtxt('data_pulse1_b0.02_d30.txt')
time_flux1, P1_1, P2_1, Z_1, PO4_1, Psupply1 = data1.T
axs[0,0].plot(time_flux1, Psupply1, label=r'$P_{supply}$', color="darkgray")
axs[0,0].plot(time_flux1, P1_1, label=r'$P_1$', color="chartreuse")
axs[0,0].plot(time_flux1, P2_1, label=r'$P_2$', color="green")
axs[0,0].plot(time_flux1, Z_1, label=r'$Z$', color="aqua")
axs[0,0].plot(time_flux1, PO4_1, label=r'$PO4$', color="magenta")
axs[0,0].set_xlabel('Time [d]', fontsize=10)
axs[0,0].set_ylabel('Masses [mmolC.m$^{-3}$]', fontsize=10)
axs[0,0].set_ylim(0, 0.6)
axs[0,0].legend(frameon=False, loc='upper left')
axs[0,0].set_title('"Low intensity b=0.02"', fontweight='bold', fontsize=12)

# Load the generated data for the second subplot
P_1 = np.loadtxt('data_senstitivitytest_b_pulse1_b0.02_d30_P1.txt')
P_2 = np.loadtxt('data_senstitivitytest_b_pulse1_b0.02_d30_P2.txt')
ratio = np.loadtxt('data_senstitivitytest_b_pulse1_b0.02_d30_ratio.txt')
P_O4 = np.loadtxt('data_senstitivitytest_b_pulse1_b0.02_d30_PO4.txt')
colors = ['#006400', '#ffffff', '#00ff00']
n_bins = 100
cmap_name = 'custom_green'
cm = plt.cm.colors.LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)
im = axs[0,1].imshow(ratio, cmap=cm, aspect='auto', origin='lower', extent=[0, (nb_time - 1) * dt, l_param[0], l_param[-1]])
axs[0,1].set_xlabel('Time [d]', fontsize=10)
axs[0,1].set_ylabel('b', fontsize=10)
axs[0,1].axhline(y=0.08, color='gray', linestyle='--')
axs[0,1].axhline(y=0.02, color='gray', linestyle='--')
axs[0,1].annotate(r'$P_2$', xy=(0.1, -0.48), xycoords='axes fraction', fontsize=12, ha='right')
axs[0,1].annotate(r'$P_1$', xy=(0.9, -0.48), xycoords='axes fraction', fontsize=12, ha='left')
cbar = plt.colorbar(im, ax=axs[0,1], orientation='horizontal', pad=0.25, shrink=0.7, extend='both', ticks=[0.4, 0.5, 0.6])
cbar.set_label(r'$R$', fontsize=10)
im.set_clim(0.4, 0.6)

# Load the generated data for the third subplot
data3 = np.loadtxt('data_pulse1_b0.08_d30.txt')
time_flux3, P1_3, P2_3, Z_3, PO4_3, Psupply3 = data3.T
axs[0,2].plot(time_flux3, Psupply3, label=r'$P_{supply}$', color="darkgray")
axs[0,2].plot(time_flux3, P1_3, label=r'$P_1$', color="chartreuse")
axs[0,2].plot(time_flux3, P2_3, label=r'$P_2$', color="green")
axs[0,2].plot(time_flux3, Z_3, label=r'$Z$', color="aqua")
axs[0,2].plot(time_flux3, PO4_3, label=r'$PO4$', color="magenta")
axs[0,2].set_xlabel('Time [d]', fontsize=10)
axs[0,2].set_ylabel('Masses [mmolC.m$^{-3}$]', fontsize=10)
axs[0,2].set_ylim(0, 0.6)
axs[0,2].set_title('"High intensity b=0.08"', fontweight='bold', fontsize=12)

# Load the generated data for the fourth subplot
Z = np.loadtxt('data_senstitivitytest_b_pulse1_b0.02_d30_Z.txt')
cm = plt.cm.colors.LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)
im = axs[1,1].imshow(Z, cmap='Greys', aspect='auto', origin='lower', extent=[0, (nb_time - 1) * dt, l_param[0], l_param[-1]])
axs[1,1].set_xlabel('Time [d]', fontsize=10)
axs[1,1].set_ylabel('b', fontsize=10)
axs[1,1].axhline(y=0.08, color='gray', linestyle='--')
axs[1,1].axhline(y=0.02, color='gray', linestyle='--')
cbar = plt.colorbar(im, ax=axs[1,1], orientation='horizontal', pad=0.25, shrink=0.7)
cbar.set_label(r'$Z$', fontsize=10)
im.set_clim(0.4, 0.5)

# Adjust layout
plt.tight_layout()
plt.show()

axs[1, 2].axis('off')
axs[1, 0].axis('off')

axs[0, 0].set_position([axs[0, 0].get_position().x0, axs[0, 0].get_position().y0 - 0.15, axs[0, 0].get_position().width, axs[0, 0].get_position().height])
axs[0, 2].set_position([axs[0, 2].get_position().x0, axs[0, 2].get_position().y0 - 0.15, axs[0, 2].get_position().width, axs[0, 2].get_position().height])


plt.savefig('../outputs/dynamicaltest1999d.pdf', format='pdf')

####### 90 days simulations ########

# Charger les données et définir les paramètres communs
dt = 0.2
end_time_flux = 90
time_flux2 = np.arange(0, end_time_flux, dt)
nb_time = len(time_flux2)
min_param = 0.01
max_param = 0.1
l_param = np.linspace(min_param, max_param, 10)
colors = ['#006400', '#ffffff', '#00ff00']
n_bins = 100
cmap_name = 'custom_green'
cm = plt.cm.colors.LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)
plt.rc('font', size=7)

# Créer la grille de sous-graphiques
fig = plt.figure(figsize=(10, 7))
gs = GridSpec(3, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1, 1]) 

# Charger les données pour les graphiques de la première rangée
data1 = np.loadtxt('data_pulse1_b0.08_d90.txt')
time_flux1, P1_1, P2_1, Z_1, PO4_1, Psupply1 = data1.T

data2 = np.loadtxt('data_pulse2_b0.08_d90.txt')
time_flux2, P1_2, P2_2, Z_2, PO4_2, Psupply2 = data2.T

data3 = np.loadtxt('data_pulse3_b0.08_d90.txt')
time_flux3, P1_3, P2_3, Z_3, PO4_3, Psupply3 = data3.T

# Créer les sous-graphiques pour les figures de la première rangée
for i, data in enumerate([(time_flux1, P1_1, P2_1, Z_1, PO4_1, Psupply1),
                           (time_flux2, P1_2, P2_2, Z_2, PO4_2, Psupply2),
                           (time_flux3, P1_3, P2_3, Z_3, PO4_3, Psupply3)]):
    ax = plt.subplot(gs[0, i])
    time_flux, P1, P2, Z, PO4, Psupply = data
    ax.plot(time_flux, Psupply, label=r'$P_{supply}$', color="darkgray")
    ax.plot(time_flux, P1, label=r'$P_1$', color="chartreuse")
    ax.plot(time_flux, P2, label=r'$P_2$', color="green")
    ax.plot(time_flux, Z, label=r'$Z$', color="aqua")
    ax.plot(time_flux, PO4, label=r'$PO4$', color="magenta")
    ax.set_xlabel('Time [d]', fontsize=10)
    ax.set_ylabel('Masses [mmolCm$^{-3}$]', fontsize=10)
    ax.set_ylim(0, 0.65)
    ax.set_xlim(0, 90)  # Définir les limites d'axe x
    if i == 0:
        ax.legend(frameon=True, loc='upper left')  # Légende uniquement pour la première figure
        ax.set_title('"One pulse"', fontweight='bold', fontsize=12)
    if i == 1:
        ax.set_title('"Two pulses"', fontweight='bold', fontsize=12)
    if i == 2:
        ax.set_title('"Three pulses"', fontweight='bold', fontsize=12)

# Créer les sous-graphiques pour les figures de la deuxième rangée
for i in range(3):
    ax = plt.subplot(gs[1, i])
    ratio = np.loadtxt(f'data_senstitivitytest_b_pulse{i+1}_b0.08_d90_ratio.txt')
    im = ax.imshow(ratio, cmap=cm, aspect='auto', origin='lower', extent=[0, (nb_time-1)*dt, l_param[0], l_param[-1]])
    ax.set_xlabel('Time [d]', fontsize=10)
    ax.set_ylabel('b [mmolCm$^{-3}d^{-1}$]', fontsize=10)
    ax.axhline(y=0.08, color='gray', linestyle='--')
    ax.annotate(r'$P_2$', xy=(0.1, -0.49), xycoords='axes fraction', fontsize=12, ha='right')
    ax.annotate(r'$P_1$', xy=(0.9, -0.49), xycoords='axes fraction', fontsize=12, ha='left')
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.25, shrink=0.7, extend='both', ticks=[0.4, 0.5, 0.6])
    cbar.set_label(r'$R$ []', fontsize=10)
    im.set_clim(0.4, 0.6)  # Utiliser set_clim() sur l'objet AxesImage
    ax.set_xlim(0, 90)  # Définir les limites d'axe x
    
# Créer les sous-graphiques pour les figures de la troisieme rangée
for i in range(3):
    ax = plt.subplot(gs[2, i])
    ratio = np.loadtxt(f'data_senstitivitytest_b_pulse{i+1}_b0.08_d90_Z.txt')
    im = ax.imshow(ratio, cmap='Greys', aspect='auto', origin='lower', extent=[0, (nb_time-1)*dt, l_param[0], l_param[-1]])
    ax.set_xlabel('Time [d]', fontsize=10)
    ax.set_ylabel('b [mmolCm$^{-3}d^{-1}$]', fontsize=10)
    ax.axhline(y=0.08, color='gray', linestyle='--')
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.25, shrink=0.7)
    cbar.set_label(r'$Z$ [mmolCm$^{-3}$]', fontsize=10)
    im.set_clim(0.4, 0.6)  # Utiliser set_clim() sur l'objet AxesImage
    ax.set_xlim(0, 90)  # Définir les limites d'axe x
    
plt.tight_layout()
plt.show()

plt.savefig('../outputs/dynamicaltest90d.pdf', format='pdf')