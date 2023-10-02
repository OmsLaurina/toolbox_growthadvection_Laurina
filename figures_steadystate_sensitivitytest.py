#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 11:52:08 2023

@author: loms
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker
from matplotlib.cm import ScalarMappable
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scienceplots
from f_monod_hollingII import f_monod
plt.style.use(['science','no-latex'])
plt.close('all')

#### SCENARIOS #####

# Set up
name_pdf = "../outputs/figure_scenarios.pdf"
name_param = r'$P_{\mathrm{supply}}$ [mmolCm$^{-3}d^{-1}$]'
test = 0.01

n= 100
min_param = 0.01
max_param = 0.1
l_param = np.linspace(min_param, max_param, n)

# Load data and create plots for each scenario
# scenarios = ['nograzing', 'equalgrazing', 'diffgrazing']
# titles = ['A. "No grazing"', 'B. "Equal grazing"', 'C. "Differential grazing"'] 

scenarios = ['diffgrazing']
titles = ["Differential grazing"]

# Create subplots with GridSpec
fig = plt.figure(figsize=(8, 6))
# gs = GridSpec(2, 3, height_ratios=[1, 1], width_ratios=[1, 1, 1], hspace=0.5, wspace=0.7)
gs = GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1], hspace=0.5, wspace=0.7)

# Upper row of subplots
for i, scenario in enumerate(scenarios):
    P_1 = np.loadtxt(f'data_senstitivitytest_Psupp_{scenario}_P1.txt')
    P_2 = np.loadtxt(f'data_senstitivitytest_Psupp_{scenario}_P2.txt')
    ratio = np.loadtxt(f'data_senstitivitytest_Psupp_{scenario}_ratio.txt')
    P_O4 = np.loadtxt(f'data_senstitivitytest_Psupp_{scenario}_PO4.txt')

    ax = fig.add_subplot(gs[0, i])
    if i == 0:
        sc = ax.scatter(l_param, ratio, c=P_1+P_2, cmap='Wistia', vmin=0, vmax=100)
    else : 
        sc = ax.scatter(l_param, ratio, c=P_1+P_2, cmap='Wistia', vmin=0.22, vmax=0.45)
        
    ax.plot(l_param, ratio, color='black')
    ax.axvline(x=0.01, color='gray', linestyle='--')
    if scenario == 'diffgrazing':
        ax.axvline(x=0.05, color='red', linestyle='--')
    ax.grid()
    ax.set_xlabel(name_param, fontsize=10)
    ax.set_ylabel(r'$R$', fontsize=10)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_title(titles[i], fontweight='bold', fontsize=12)

    cbar_ax = ax.inset_axes([1.05, 0.2, 0.05, 0.6])
    format_str = '%.0f' if i == 0 else '%.2f'
    cbar = fig.colorbar(sc, cax=cbar_ax, format=format_str,ticks=MaxNLocator(nbins=3)) 
    cbar.set_label(r'$P_1+P_2$', labelpad=2)
    

# Lower row of subplots
for i, scenario in enumerate(scenarios):
    data = np.loadtxt(f'../outputs/data_{scenario}_{test}.txt')
    time, P1, P2, Z, PO4 = data.T
    ax = fig.add_subplot(gs[1, i])
    ax.plot(time, P1, label=r'$P_1$', color="chartreuse")
    ax.plot(time, P2, label=r'$P_2$', color="green")
    ax.plot(time, Z, label=r'$Z$', color="aqua")
    ax.plot(time, PO4, label=r'$PO4$', color="magenta")
    ax.set_xlabel('Time [d]', fontsize=10)
    ax.set_ylabel('Masses [mmolC.m$^{-3}$]', fontsize=10)
    ax.set_xticks([0, 500, 1000, 1500, 2000])
    y_limit = (10 if i == 0 else 0.9)
    ax.set_ylim(0, y_limit)
    
    if i == 0:
        legend = ax.legend(frameon=False, loc='upper left') 
        for label in legend.get_texts():
            label.set_fontsize(7)

# Adjust layout and save the figure
plt.tight_layout()
plt.savefig(name_pdf, format='pdf')

########################## 

name_pdf = "../outputs/figure_diffgrazing.pdf"
# Create a figure with GridSpec
fig = plt.figure(figsize=(12,3))
gs = GridSpec(1, 3, width_ratios=[1, 1, 1], wspace=0.3)

# Upper row of subplots
P_1 = np.loadtxt(f'data_senstitivitytest_Psupp_{scenario}_P1.txt')
P_2 = np.loadtxt(f'data_senstitivitytest_Psupp_{scenario}_P2.txt')
ratio = np.loadtxt(f'data_senstitivitytest_Psupp_{scenario}_ratio.txt')
P_O4 = np.loadtxt(f'data_senstitivitytest_Psupp_{scenario}_PO4.txt')

ax = fig.add_subplot(gs[1])
sc = ax.scatter(l_param, ratio, c=P_1 + P_2, cmap='Wistia')

ax.plot(l_param, ratio, color='black')
ax.axvline(x=0.01, color='gray', linestyle='--',linewidth=2)
ax.axvline(x=0.1, color='gray', linestyle='--',linewidth=2)
if scenario == 'diffgrazing':
    ax.axvline(x=0.05, color='red', linestyle='--')
ax.grid()
ax.set_xlabel(name_param, fontsize=10)
ax.set_ylabel('R []', fontsize=10)
ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
# ax.set_title(titles[i], fontweight='bold', fontsize=12)

# cbar_ax = ax.inset_axes([1.05, 0.2, 0.05, 0.6])
# cbar = fig.colorbar(sc, cax=cbar_ax, format=format_str, ticks=MaxNLocator(nbins=3))
# cbar.set_label(r'$P_1+P_2$ [mmolCm$^{-3}$]', labelpad=2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("top", size="5%", pad=0.2)  # Ajustez la taille et le pad selon vos besoins
cbar = fig.colorbar(sc, cax=cax, orientation="horizontal", format=format_str, ticks=MaxNLocator(nbins=3))
cbar.ax.xaxis.set_label_position('top')
cbar.set_label(r'$P_1+P_2$ [mmolCm$^{-3}$]', labelpad=2)

# Lower row of subplots
data = np.loadtxt(f'../outputs/data_{scenario}_0.01.txt')
time, P1, P2, Z, PO4 = data.T
ax = fig.add_subplot(gs[0])
ax.plot(time, P1, label=r'$P_1$', color="chartreuse")
ax.plot(time, P2, label=r'$P_2$', color="green")
ax.plot(time, Z, label=r'$Z$', color="aqua")
ax.plot(time, PO4, label=r'$PO4$', color="magenta")
ax.set_xlabel('Time [d]', fontsize=10)
ax.set_ylabel(r'Masses [mmolCm$^{-3}$]', fontsize=10)
ax.set_xticks([0, 500, 1000])
ax.set_xlim(0,1000)

legend = ax.legend(frameon=False, loc='upper right')
for label in legend.get_texts():
    label.set_fontsize(7)
    
data = np.loadtxt(f'../outputs/data_{scenario}_0.1.txt')
time, P1, P2, Z, PO4 = data.T
ax = fig.add_subplot(gs[2])
ax.plot(time, P1, label=r'$P_1$', color="chartreuse")
ax.plot(time, P2, label=r'$P_2$', color="green")
ax.plot(time, Z, label=r'$Z$', color="aqua")
ax.plot(time, PO4, label=r'$PO4$', color="magenta")
ax.set_xlabel('Time [d]', fontsize=10)
ax.set_ylabel(r'Masses [mmolCm$^{-3}$]', fontsize=10)
ax.set_xticks([0, 500, 1000])
ax.set_xlim(0,1000)
    
plt.tight_layout()
plt.savefig(name_pdf, format='pdf')
plt.show()

##### MONOD ######

plt.rc('font', size=7)

plt.figure(3)

name_pdf_monod = '../outputs/monod.pdf'
test = 0.1
data = np.loadtxt(f'../outputs/data_{test}.txt')
time, P1, P2, Z, PO4 = data.T

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

N_theo_PO4 = np.linspace(0, 6, len(time))
f_monod(N_theo_PO4, arg['umax1'], arg['kP1'],"chartreuse")
f_monod(N_theo_PO4, arg['umax2'], arg['kP2'],"green")
# f_monod(P_O4, arg['umax1'], arg['kP1'],"gray")
# f_monod(P_O4, arg['umax2'], arg['kP2'],"lightgray")
plt.xlabel('PO4 [mmolCm\u207B\u00B3]')
plt.ylabel(r'$\mu\ [d^{-1}]$')
plt.ylim(0,2)
plt.axvline(x=1.378, color='red', linestyle='--', label="[PO4] observed")
plt.axvline(x=min(P_O4), color='lightcoral', linestyle='--', label="Min PO4")
plt.axvline(x=max(P_O4), color='lightcoral', linestyle='--', label="Max PO4")
plt.fill_betweenx(np.linspace(0, max(plt.ylim()), n), min(P_O4), max(P_O4), color='red', alpha=0.3, label='Range')

plt.legend(handles=[plt.Line2D([0], [0], color='chartreuse', label=r'$P_1$'),
                plt.Line2D([0], [0], color='green', label=r'$P_2$'),
                plt.Line2D([0], [0], color='red', linestyle='--', label='Observed [PO4]'),
                plt.Line2D([0], [0], color='lightcoral', linestyle='--', label='Modeled [PO4] range')
                ])

plt.savefig(name_pdf_monod, format='pdf')


# #### PROPORTIONS ####

# # Load data from text files
# with open("prop_p1_param1_inf_0.05.txt", "r") as f_prop_p1_inf_005:
#     prop_p1_inf_005 = [float(line.strip().split('\t')[0]) for line in f_prop_p1_inf_005.readlines()]

# with open("prop_p2_param1_inf_0.05.txt", "r") as f_prop_p2_inf_005:
#     prop_p2_inf_005 = [float(line.strip().split('\t')[0]) for line in f_prop_p2_inf_005.readlines()]

# with open("prop_p1_param1_sup_0.05.txt", "r") as f_prop_p1_sup_005:
#     prop_p1_sup_005 = [float(line.strip().split('\t')[0]) for line in f_prop_p1_sup_005.readlines()]

# with open("prop_p2_param1_sup_0.05.txt", "r") as f_prop_p2_sup_005:
#     prop_p2_sup_005 = [float(line.strip().split('\t')[0]) for line in f_prop_p2_sup_005.readlines()]

# # Extract param1 values from the same files
# with open("prop_p1_param1_inf_0.05.txt", "r") as f_prop_p1_inf_005:
#     param1_inf_005 = [float(line.strip().split('\t')[1]) for line in f_prop_p1_inf_005.readlines()]

# with open("prop_p1_param1_sup_0.05.txt", "r") as f_prop_p1_sup_005:
#     param1_sup_005 = [float(line.strip().split('\t')[1]) for line in f_prop_p1_sup_005.readlines()]

# fig = plt.figure(figsize=(6, 2))

# ax1 = fig.add_subplot(1, 2, 1)
# ax2 = fig.add_subplot(1, 2, 2)

# # Plot for param1 < 0.05

# # Plot for param1 < 0.05
# indices_inf_005 = np.arange(len(prop_p1_inf_005))
# bar_width = 0.4

# ax1.bar(indices_inf_005 - bar_width / 2, prop_p1_inf_005, width=bar_width, label=r'$P_1$', color='chartreuse', alpha=0.7)
# ax1.bar(indices_inf_005 + bar_width / 2, prop_p2_inf_005, width=bar_width, label=r'$P_2$', color='green', alpha=0.7)

# ax1.set_xticks(indices_inf_005)
# ax1.set_xticklabels([round(val, 3) for val in param1_inf_005], rotation=90)  # Set param1 values as x-axis labels
# ax1.set_xlabel(r'$P_{supply}$')
# ax1.set_ylabel('Proportion (%)')
# ax1.legend()
# ax1.set_title('Low-nutrient water')
# ax1.set_xticks(indices_inf_005[::5])

# # Plot for param1 >= 0.05
# indices_sup_005 = np.arange(len(prop_p1_sup_005))

# ax2.bar(indices_sup_005 - bar_width / 2, prop_p1_sup_005, width=bar_width, label='p1', color='chartreuse', alpha=0.7)
# ax2.bar(indices_sup_005 + bar_width / 2, prop_p2_sup_005, width=bar_width, label='p2', color='green', alpha=0.7)

# # param1_sup_005 = [round(float(line.strip().split('\t')[1]), 2) for line in f_prop_p1_sup_005.readlines()]
# ax2.set_xticks(indices_sup_005)
# ax2.set_xticklabels([round(val, 3) for val in param1_sup_005], rotation=90)  # Set param1 values as x-axis labels
# ax2.set_xlabel(r'$P_{supply}$')
# ax2.set_title('High-nutrient water')
# ax2.set_xticks(indices_sup_005[::5])

# # Remove ticks from the y-axis of the second plot
# ax2.set_yticks([])
# plt.tight_layout()
# plt.show()
# plt.savefig('percentage.pdf', format='pdf')