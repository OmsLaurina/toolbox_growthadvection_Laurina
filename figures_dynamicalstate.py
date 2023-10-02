#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 17:06:26 2023

@author: loms
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


plt.figure(1)
# Load the generated data
data = np.loadtxt('data_pulse3_b0.08_d90.txt')
time_flux, P1, P2, Z, PO4, Psupply = data.T

name_pdf = 'data_pulse3_b0.08_d90.pdf'
plt.rc('font', size=7) 
plt.plot(time_flux,Psupply, label=r'$P_{supply}$', color="darkgray")
plt.plot(time_flux, P1, label=r'$P_1$', color="chartreuse")
plt.plot(time_flux, P2, label=r'$P_2$', color="green")
plt.plot(time_flux, Z, label= r'$Z$', color="aqua")
plt.plot(time_flux,PO4, label=r'$PO4$', color="magenta")
plt.xlabel('Time [day]',fontsize=10)
plt.ylabel('Masses  [mmolC.m\u207B\u00B3]',fontsize=10)
plt.ylim(0,0.6)
plt.legend(frameon=True)
legend = plt.legend(frameon=True, bbox_to_anchor=(1.05, 1), loc='upper left')
# plt.text(0.15, 0.9, f'Pmoy = {moy_flux}',transform=plt.gcf().transFigure, ha='center', color='red')
# plt.text(0.5, 0.9, f'b = {b}',transform=plt.gcf().transFigure, ha='center', color='red')
# plt.text(0.8, 0.9, f'T = {T*dt}',transform=plt.gcf().transFigure, ha='center', color='red')
plt.title(' "Three pulses" ',fontweight='bold',fontsize=12)

plt.figure(2)
# Load the generated data
data = np.loadtxt('data_pulse2_b0.08_d90.txt')
time_flux, P1, P2, Z, PO4, Psupply = data.T
name_pdf = 'data_pulse2_b0.08_d90.pdf'
plt.rc('font', size=7) 
plt.plot(time_flux,Psupply, label=r'$P_{supply}$', color="darkgray")
plt.plot(time_flux, P1, label=r'$P_1$', color="chartreuse")
plt.plot(time_flux, P2, label=r'$P_2$', color="green")
plt.plot(time_flux, Z, label= r'$Z$', color="aqua")
plt.plot(time_flux,PO4, label=r'$PO4$', color="magenta")
plt.xlabel('Time [day]',fontsize=10)
plt.ylabel('Masses  [mmolC.m\u207B\u00B3]',fontsize=10)
plt.ylim(0,0.6)
plt.legend(frameon=True)
legend = plt.legend(frameon=True, bbox_to_anchor=(1.05, 1), loc='upper left')
# plt.text(0.15, 0.9, f'Pmoy = {moy_flux}',transform=plt.gcf().transFigure, ha='center', color='red')
# plt.text(0.5, 0.9, f'b = {b}',transform=plt.gcf().transFigure, ha='center', color='red')
# plt.text(0.8, 0.9, f'T = {T*dt}',transform=plt.gcf().transFigure, ha='center', color='red')
plt.title(' "Two pulses" ',fontweight='bold',fontsize=12)

plt.figure(3)
# Load the generated data
data = np.loadtxt('data_pulse1_b0.08_d90.txt')
time_flux, P1, P2, Z, PO4, Psupply = data.T
name_pdf = 'data_pulse1_b0.08_d90.pdf'
plt.rc('font', size=7) 
plt.plot(time_flux,Psupply, label=r'$P_{supply}$', color="darkgray")
plt.plot(time_flux, P1, label=r'$P_1$', color="chartreuse")
plt.plot(time_flux, P2, label=r'$P_2$', color="green")
plt.plot(time_flux, Z, label= r'$Z$', color="aqua")
plt.plot(time_flux,PO4, label=r'$PO4$', color="magenta")
plt.xlabel('Time [day]',fontsize=10)
plt.ylabel('Masses  [mmolC.m\u207B\u00B3]',fontsize=10)
plt.ylim(0,0.6)
plt.legend(frameon=True)
legend = plt.legend(frameon=True, bbox_to_anchor=(1.05, 1), loc='upper left')
# plt.text(0.15, 0.9, f'Pmoy = {moy_flux}',transform=plt.gcf().transFigure, ha='center', color='red')
# plt.text(0.5, 0.9, f'b = {b}',transform=plt.gcf().transFigure, ha='center', color='red')
# plt.text(0.8, 0.9, f'T = {T*dt}',transform=plt.gcf().transFigure, ha='center', color='red')
plt.title(' "One pulse" ',fontweight='bold',fontsize=12)


plt.figure(4)
# Load the generated data
data = np.loadtxt('data_pulse1_b0.08_d30.txt')
time_flux, P1, P2, Z, PO4, Psupply = data.T
name_pdf = 'data_pulse1_b0.08_d30.pdf'
plt.rc('font', size=7) 
plt.plot(time_flux,Psupply, label=r'$P_{supply}$', color="darkgray")
plt.plot(time_flux, P1, label=r'$P_1$', color="chartreuse")
plt.plot(time_flux, P2, label=r'$P_2$', color="green")
plt.plot(time_flux, Z, label= r'$Z$', color="aqua")
plt.plot(time_flux,PO4, label=r'$PO4$', color="magenta")
plt.xlabel('Time [day]',fontsize=10)
plt.ylabel('Masses  [mmolC.m\u207B\u00B3]',fontsize=10)
plt.ylim(0,0.6)
plt.legend(frameon=True)
legend = plt.legend(frameon=True, bbox_to_anchor=(1.05, 1), loc='upper left')
# plt.text(0.15, 0.9, f'Pmoy = {moy_flux}',transform=plt.gcf().transFigure, ha='center', color='red')
# plt.text(0.5, 0.9, f'b = {b}',transform=plt.gcf().transFigure, ha='center', color='red')
# plt.text(0.8, 0.9, f'T = {T*dt}',transform=plt.gcf().transFigure, ha='center', color='red')
plt.title(' "High intensity b=0.08" ',fontweight='bold',fontsize=12)


plt.figure(5)
# Load the generated data
data = np.loadtxt('data_pulse1_b0.02_d30.txt')
time_flux, P1, P2, Z, PO4, Psupply = data.T
name_pdf = 'data_pulse1_b0.02_d30.pdf'
plt.rc('font', size=7) 
plt.plot(time_flux,Psupply, label=r'$P_{supply}$', color="darkgray")
plt.plot(time_flux, P1, label=r'$P_1$', color="chartreuse")
plt.plot(time_flux, P2, label=r'$P_2$', color="green")
plt.plot(time_flux, Z, label= r'$Z$', color="aqua")
plt.plot(time_flux,PO4, label=r'$PO4$', color="magenta")
plt.xlabel('Time [day]',fontsize=10)
plt.ylabel('Masses  [mmolC.m\u207B\u00B3]',fontsize=10)
plt.ylim(0,0.6)
plt.legend(frameon=True)
legend = plt.legend(frameon=True, bbox_to_anchor=(1.05, 1), loc='upper left')
# plt.text(0.15, 0.9, f'Pmoy = {moy_flux}',transform=plt.gcf().transFigure, ha='center', color='red')
# plt.text(0.5, 0.9, f'b = {b}',transform=plt.gcf().transFigure, ha='center', color='red')
# plt.text(0.8, 0.9, f'T = {T*dt}',transform=plt.gcf().transFigure, ha='center', color='red')
plt.title(' "Low intensity b=0.02" ',fontweight='bold',fontsize=12)