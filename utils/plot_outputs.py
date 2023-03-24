"""
Plot the outputs of the growth function : 
    1. Temporal evolution of variables
    2. Monod function
    3. Portrait de phase
"""
import matplotlib.pyplot as plt
import numpy as npy
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import scienceplots
import sys
sys.path.append('../')

plt.style.use(['science','no-latex'])
plt.close('all')

""" Set up """

# Timestep (same than timestep of the growth function : by default is 0.2)
dt = 0.2
# Choose the time at which the simulation end (days)
end_time = 2000
# create a time vector from 0 to end_time, seperated by a timestep
time = npy.arange(0,end_time,dt)

# Externe input of nutrients 
Psupply_moy = 0.1 #0.05 = threshold where P2 became dominant
#Psupply constant
Psupply =[Psupply_moy] * len(time)

# Call the function
from growth_model_2P1Z_v10 import growth_model_2P1Z_v10
[P1,P2,Z,PO4,Export,u1,u2,g1,g2,arg]=growth_model_2P1Z_v10(Psupply, time)

Export.append(Export[-1])

""" Figures """

from pylab import rcParams
rcParams['figure.figsize'] = 4,3.5

from matplotlib.backends.backend_pdf import PdfPages
with PdfPages('../outputs/model_V10/outputs_analysis_v10.pdf') as pdf:

    # Figure 1 : temporal evolution of biomasses
    plt.figure(1)
    plt.plot(time, P1, label='P1', color="chartreuse")
    plt.plot(time, P2, label='P2', color="green")
    plt.plot(time, Z, label='Zoo', color="aqua")
    plt.plot(time,PO4, label='PO4', color="magenta")
    plt.xlabel('Time')
    plt.ylabel('mmolC.m^-3')
    plt.legend(frameon=True) 
    plt.title('model outputs (concentration over time)')
    pdf.savefig()
    
    #Budget
    plt.figure(2)
    plt.plot(time, Psupply, label='Psupply', color="chartreuse")
    plt.plot(time, Export, label='Export', color="green")
    plt.xlabel('Time')
    plt.ylabel('mmolC.m^-3.d^-1')
    plt.legend(frameon=True) 
    plt.title('Budget')
    pdf.savefig()
    
    #Figure 2 : Monod function and Holling type II response
    from utils.f_monod_hollingII import f_monod, f_hollingII
    
    ## MONOD
    plt.figure(3)
    N_theo_PO4 = np.linspace(0, 10, len(time))
    f_monod(N_theo_PO4, arg['umax1'], arg['kP1'],"chartreuse")
    f_monod(N_theo_PO4, arg['umax2'], arg['kP2'],"green")
    f_monod(PO4, arg['umax1'], arg['kP1'],"blue")
    f_monod(PO4, arg['umax2'], arg['kP2'],"red")
    
    plt.xlabel('[PO4] mmolC.m³')
    plt.ylabel('\u03BC d^{-1}')
    plt.legend(['P_1_theo', 'P_2_theo', 'P_1_mod', 'P_2_mod'])
    pdf.savefig()
      
    ## HOLLINGII
    plt.figure(4)
    P1_theo = np.linspace(0, 5, len(time))
    P2_theo = np.linspace(0, 5, len(time))
    SUM_P = P1_theo + P2_theo
    
    P1_array = np.array(P1)
    P2_array = np.array(P2)
    SUM_P_array = P1_array+P2_array
    
    f_hollingII(P1_theo, SUM_P, arg['gmax1'], arg['kZ1'],"chartreuse")
    f_hollingII(P2_theo, SUM_P, arg['gmax2'], arg['kZ2'],"green")
    f_hollingII(P1, SUM_P_array, arg['gmax1'], arg['kZ1'],"blue")
    f_hollingII(P2, SUM_P_array, arg['gmax2'], arg['kZ2'],"red")
    
    plt.xlabel('Phytoplancton mmolC.m³')
    plt.ylabel('g d^{-1}')
    plt.legend(['Z-->P1_theo', 'Z-->P2_theo', 'Z-->P1_mod', 'Z-->P2_mod'])
    pdf.savefig()
    
    # # Figure 3 : Portrait de phase
    # # 1. P1
    # plt.figure(5)
    # ax1 = plt.axes(projection ='3d')
    
    # # Calcul des dérivées de P1 et PO4 par rapport au temps
    # dp1_dt = np.gradient(P1, dt)
    # dZ_dt = np.gradient(Z, dt)
    # dpo4_dt = np.gradient(PO4, dt)
    
    # # Tracé du portrait de phase
    # ax1.plot(P1,Z,PO4,label='model trajectory')
    # ax1.set_xlabel('P1')
    # ax1.set_ylabel('Z')
    # ax1.set_zlabel('PO4')
    # ax1.legend()
    
    # # Tracé du champ de vecteur sur le portrait de phase
    # ax1.quiver(P1, Z, PO4, dp1_dt, dZ_dt, dpo4_dt, edgecolor='red',length=0.5, arrow_length_ratio=0.3, linewidth=0.5)
    # ax1.view_init(23, -151)
    # pdf.savefig()
    
    # # 2. P2
    # plt.figure(6)
    # ax2 = plt.axes(projection ='3d')
    
    # dp2_dt = np.gradient(P2, dt)
    # dZ_dt = np.gradient(Z, dt)
    # dpo4_dt = np.gradient(PO4, dt)
    
    # ax2.plot(P2,Z,PO4,label='model trajectory')
    # ax2.set_xlabel('P2')
    # ax2.set_ylabel('Z')
    # ax2.set_zlabel('PO4')
    # ax2.legend()
    
    # ax2.quiver(P2, Z, PO4, dp2_dt, dZ_dt, dpo4_dt, edgecolor='red',length=0.5, arrow_length_ratio=0.3, linewidth=0.5)
    # ax2.view_init(23, -151)
    # pdf.savefig()