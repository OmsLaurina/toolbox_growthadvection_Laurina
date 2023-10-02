#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 11:22:10 2023

@author: loms
"""

import matplotlib.pyplot as plt
import numpy as np
import scienceplots
from f_monod_hollingII import f_monod


test = '0.1'

# Load data from the saved file
data = np.loadtxt(f'../outputs/data_{test}.txt')
time, P1, P2, Z, PO4 = data.T

plt.style.use(['science', 'no-latex'])
plt.close('all')

# Figure 1 : Temporal evolution of biomasses
plt.figure(1)
plt.rc('font', size=7)
plt.plot(time, P1, label=r'$P_1$', color="chartreuse")
plt.plot(time, P2, label=r'$P_2$', color="green")
plt.plot(time, Z, label=r'$Z$', color="aqua")
plt.plot(time, PO4, label=r'$PO4$', color="magenta")

plt.xlabel('Time [day]', fontsize=10)
plt.ylabel('Masses [mmolC.m\u207B\u00B3]', fontsize=10)
plt.xticks([0, 500, 1000, 1500, 2000])
legend = plt.legend(frameon=True, loc='upper right')
plt.gca().add_artist(legend)
figure_filename = f'../outputs/tempevol_{test}.pdf'
plt.savefig(figure_filename, format='pdf') 


# Figure 2 :  Monod curves

arg = {
    'umax1': 1.9872,
    'umax2': 2.7648,
    'kP1': 1,
    'kP2': 3,
    'gmax1':3.89,
    'gmax2':0.43,
    'kZ1': 5,
    'kZ2': 20,
}

plt.figure(2)
N_theo_PO4 = np.linspace(0, 6, len(time))
f_monod(N_theo_PO4, arg['umax1'], arg['kP1'], "chartreuse")
f_monod(N_theo_PO4, arg['umax2'], arg['kP2'], "green")
f_monod(PO4, arg['umax1'], arg['kP1'], "gray")
f_monod(PO4, arg['umax2'], arg['kP2'], "lightgray")
plt.xlabel('PO4 [mmolC.m\u207B\u00B3]')
plt.ylabel(r'$\mu [d^{-1}]$')
plt.axvline(x=1.378, color='red', linestyle='--', label="[PO4] standard")
plt.axvline(x=PO4[-1], color='lightcoral', linestyle='--', label="[PO4] modelized")
plt.legend(handles=[plt.Line2D([0], [0], color='chartreuse', label=r'$P_1$'),
                    plt.Line2D([0], [0], color='green', label=r'$P_2$'),
                    plt.Line2D([0], [0], color='red', linestyle='--', label='[PO4] standard'),
                    plt.Line2D([0], [0], color='lightcoral', linestyle='--', label='[PO4] modelized')

                ])
figure_filename = f'../outputs/monod_{test}.pdf'
plt.savefig(figure_filename, format='pdf')


# Figure 3 : Histogramme

# Charger les données pour Psupply = 0.01
data_001 = np.loadtxt('../outputs/data_0.01.txt')
time_001, P1_001, P2_001, Z_001, PO4_001 = data_001.T

Tot_001 = P1_001[-1] + P2_001[-1]
prop_P1_001 = (P1_001[-1] / Tot_001) * 100
prop_P2_001 = (P2_001[-1] / Tot_001) * 100

# Charger les données pour Psupply = 0.1
data_01 = np.loadtxt('../outputs/data_0.1.txt')
time_01, P1_01, P2_01, Z_01, PO4_01 = data_01.T

Tot_01 = P1_01[-1] + P2_01[-1]
prop_P1_01 = (P1_01[-1] / Tot_01) * 100
prop_P2_01 = (P2_01[-1] / Tot_01) * 100

# Couleurs pour les barres P1 et P2
color_P1 = 'Chartreuse'
color_P2 = 'Green'

# Créer le graphique
fig, ax = plt.subplots(figsize=(10, 6))

# Créer les barres pour Psupply = 0.01
ax.bar(0, prop_P1_001, width=0.4, label='prop_P1 (Psupply = 0.01)', color=color_P1)
ax.bar(0.5, prop_P2_001, width=0.4, label='prop_P2 (Psupply = 0.01)', color=color_P2)

# Créer les barres pour Psupply = 0.1
ax.bar(1.5, prop_P1_01, width=0.4, label='prop_P1 (Psupply = 0.1)', color=color_P1)
ax.bar(2, prop_P2_01, width=0.4, label='prop_P2 (Psupply = 0.1)', color=color_P2)

# Ajouter des étiquettes et une ligne verticale de séparation
ax.set_xticks([0.25, 1.75])
ax.set_xticklabels(['Psupply = 0.01', 'Psupply = 0.1'])
ax.axvline(x=1, color='black', linestyle='--')
ax.set_ylabel('Proportions (%)')
ax.legend()
plt.show()

figure_filename = '../outputs/histogramme.pdf'
plt.savefig(figure_filename, format='pdf')

