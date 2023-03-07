#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical analysis of the growth model

    - define differential function
    - solve them in time
    - plot the result
     
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import scienceplots

plt.style.use(['science','no-latex'])
plt.close('all')

# Set up

end_time = 600
dt = 0.2
time=list(np.arange(0,end_time,dt))
p10=0.6;p20=0.1;z0=0.6;po40=0.5

umax1=1.9872  # maximum growth rates of P1 [d**{-1}]
umax2=2.7648  # maximum growth rates of P2 [d**{-1}]
gmax1=0.5     # maximum grazing rates of Z on P1 [d**{-1}]
gmax2=0.4     # maximum grazing rates of Z on P2 [d**{-1}]
kP1=1         # half-saturation constant for P1 on PO4 [mmolC m**{-3}]
kP2=2	      # half-saturation constant for P2 on PO4 [mmolC m**{-3}]
kZ=5		  # half-saturation constant for Z on P [mmolC m**{-3}]
mP=0		  # P2 mortality rate (default 0 ie no P2 sinking) [d**{-1}]
mZ=0.07 	  # Z2 quadratic mortality rate [mmolC**{-1} m**{3} d**{-1}]
eZ=0.1        # zoo excretion rate (Z1 and Z2) [d**{-1}]
gamma1=0.7   # conversion factor from P1 to Z [/]
gamma2=0.7    # conversion factor from P2 to Z [/]
epsilon=0.25  # fraction of Z excretion that is available as regenerated PO4 [/]
Psupply = 0.01


# Define the function of growth model

def growth_model(t,u):

    global umax1, umax2, gmax1, gmax2, kP1, kP2, kZ, eZ, epsilon, mZ, gamma1, gamma2, Psupply
    
    p1=u[0];p2=u[1];z=u[2];po4=u[3] 
    
    dp1 = po4/(kP1+po4)*umax1*p1-p1/(p1+p2+kZ)*gmax1*z
    dp2 = po4/(kP2+po4)*umax2*p2-p2/(p1+p2+kZ)*gmax2*z
    dz = gamma1*p1/(p1+p2+kZ)*gmax1*z+gamma2*p2/(p1+p2+kZ)*gmax2*z-eZ*((1-gamma1)*p1/(p1+p2+kZ)*gmax1*z+(1-gamma2)*p2/(p1+p2+kZ)*gmax2*z)-mZ*z**2
    dpo4 = Psupply+epsilon*eZ*((1-gamma1)*p1/(p1+p2+kZ)*gmax1*z+(1-gamma2)*p2/(p1+p2+kZ)*gmax2*z)-po4/(kP1+po4)*umax1*p1+po4/(kP2+po4)*umax2*p2
    
    return [dp1, dp2, dz, dpo4] 

def growth_model2(u):

    global umax1, umax2, gmax1, gmax2, kP1, kP2, kZ, eZ, epsilon, mZ, gamma1, gamma2, Psupply
    
    p1=u[0];p2=u[1];z=u[2];po4=u[3] 
    
    dp1 = po4/(kP1+po4)*umax1*p1-p1/(p1+p2+kZ)*gmax1*z
    dp2 = po4/(kP2+po4)*umax2*p2-p2/(p1+p2+kZ)*gmax2*z
    dz = gamma1*p1/(p1+p2+kZ)*gmax1*z+gamma2*p2/(p1+p2+kZ)*gmax2*z-eZ*((1-gamma1)*p1/(p1+p2+kZ)*gmax1*z+(1-gamma2)*p2/(p1+p2+kZ)*gmax2*z)-mZ*z**2
    dpo4 = Psupply+epsilon*eZ*((1-gamma1)*p1/(p1+p2+kZ)*gmax1*z+(1-gamma2)*p2/(p1+p2+kZ)*gmax2*z)-po4/(kP1+po4)*umax1*p1+po4/(kP2+po4)*umax2*p2
    
    return [dp1, dp2, dz, dpo4] 

# Simplified model

def growth_model_simplify(t,u):

    global umax1, umax2, gmax1, gmax2, eZ, epsilon, mZ, gamma1, gamma2, Psupply
    
    p1=u[0];p2=u[1];z=u[2];po4=u[3] 
    
    dp1 = po4*umax1*p1-p1*gmax1*z
    dp2 = po4*umax2*p2-p2*gmax2*z
    dz = gamma1*p1*gmax1*z+gamma2*p2*gmax2*z-eZ*((1-gamma1)*p1*gmax1*z+(1-gamma2)*p2*gmax2*z)-mZ*z**2
    dpo4 = Psupply+epsilon*eZ*((1-gamma1)*p1*gmax1*z+(1-gamma2)*p2*gmax2*z)-(po4*umax1*p1+po4*umax2*p2)
    
    return [dp1, dp2, dz, dpo4] 

def growth_model2_simplify(u):

    global umax1, umax2, gmax1, gmax2, eZ, epsilon, mZ, gamma1, gamma2, Psupply
    
    p1=u[0];p2=u[1];z=u[2];po4=u[3] 
    
    dp1 = po4*umax1*p1-p1*gmax1*z
    dp2 = po4*umax2*p2-p2*gmax2*z
    dz = gamma1*p1*gmax1*z+gamma2*p2*gmax2*z-eZ*((1-gamma1)*p1*gmax1*z+(1-gamma2)*p2*gmax2*z)-mZ*z**2
    dpo4 = Psupply+epsilon*eZ*((1-gamma1)*p1*gmax1*z+(1-gamma2)*p2*gmax2*z)-(po4*umax1*p1+po4*umax2*p2)
    
    return [dp1, dp2, dz, dpo4] 

# Initial conditions
CI=[p10,p20,z0,po40]

# use function fsolve
# Cette fonction renvoit une structure complexe, qu'on note sol ici. sol est une structure qui a de nombreux 
# attributs. Parmi ces attributs, il y  notamment sol.t et sol.y, qui contiennent respectivement la liste
# des temps où les variables d'état sont données, et les variables d'état calculées pour ces temps.

sol=solve_ivp(growth_model_simplify,[0,end_time],CI,method='RK45',t_eval=time)
t=sol.t
p1=sol.y[0,:];p2=sol.y[1,:];z=sol.y[2,:];po4=sol.y[3,:]


u0=CI
Psupplymax=1;dK=0.0001;
Xeq=[];Yeq=[];
Xeq_ana=[];Yeq_ana=[]
Psupply=Psupplymax
VP=[]
KK=[]
while Psupply>dK:
    equilibrium = fsolve(growth_model2_simplify, u0) # remarque : tester CI et u0 comme valeur d'initiation de l'algo
    Xeq=Xeq+[equilibrium[0]]
    Yeq=Yeq+[equilibrium[1]]
    u0=equilibrium
    Xeq_ana=Xeq_ana+[Psupply]
    Yeq_ana=Yeq_ana+[0]
    j11 = equilibrium[2]*umax1-gmax1*equilibrium[2];
    j12 = 0;
    j13 = 0;
    j14 = 0;
    j21 = 0;
    j22 = equilibrium[3]*umax2-gmax2*equilibrium[2];
    j23 = equilibrium[3]*umax2*equilibrium[1]-gmax2*equilibrium[1];
    j24 = umax2*equilibrium[1]-equilibrium[1]*gmax2*equilibrium[2];
    j31 = gamma1*gmax1*gamma2*equilibrium[1]*gmax2*equilibrium[2]-eZ*((1-gamma1)*gmax1*equilibrium[2]+(1-gamma2)*equilibrium[1]*gmax2*equilibrium[2])-mZ*equilibrium[2]**2;
    j32 = gamma2*gmax2-eZ*(1-gamma2)*gmax2-mZ*equilibrium[2]**2;
    j33 = equilibrium[1]*gmax2*gamma2-eZ*gamma2*gmax2-2*mZ*equilibrium[2];
    j34 = 0;
    j41 = Psupply+epsilon*eZ*((1-gamma1)*gmax1*equilibrium[2]+(1-gamma2)*equilibrium[1]*gmax2*equilibrium[2])-equilibrium[3]*umax1+equilibrium[3]*equilibrium[1]*umax2;
    j42 = Psupply+epsilon*eZ*((1-gamma2)*gmax2*equilibrium[2])-equilibrium[3]*umax2;
    j43 = Psupply+epsilon*eZ*((1-gamma2)*equilibrium[1]*gmax2)-equilibrium[3]*umax2*equilibrium[1];
    j44 = Psupply+epsilon*eZ*((1-gamma2)*equilibrium[1]*gmax2*equilibrium[2])-umax2*equilibrium[1];
    J = [[j11,j12,j13,j14], [j21,j22,j23,j24], [j31,j32,j33,j34], [j41,j42,j43,j44]]
    J = np.array(J)
    vp=np.linalg.eig(J)
    VP=VP+[max(np.real(vp[0]))]
    KK=KK+[Psupply]
    Psupply=Psupply-dK


from pylab import rcParams
rcParams['figure.figsize'] = 4,3.5

from matplotlib.backends.backend_pdf import PdfPages
with PdfPages('outputs/model_V9/outputs_numerical_analysis_v9.pdf') as pdf:
    # Figures
    
    # Temporal evolution
    plt.figure(1)
    plt.plot(t, p1, label='P1', color="chartreuse")
    plt.plot(t, p2, label='P2', color="green")
    plt.plot(t, z, label='Zoo', color="aqua")
    plt.plot(t,po4, label='PO4', color="magenta")
    plt.xlabel('Time')
    plt.ylabel('mmolC.m**-3')
    plt.legend(frameon=True) 
    plt.title('model outputs (concentration over time)')
    pdf.savefig()
    
    # Portrait de phase
    
    #P1
    plt.figure(2)
    ax = plt.axes(projection ='3d')
    ax.plot(p1,po4,z)
    ax.set_xlabel('P1')
    ax.set_ylabel('PO4')
    ax.set_zlabel('Z')
    ax.view_init(23, -151)
    pdf.savefig()
    
    #P2
    plt.figure(3)
    ax2 = plt.axes(projection ='3d')
    ax2.plot(p2,po4,z)
    ax2.set_xlabel('P2')
    ax2.set_ylabel('PO4')
    ax2.set_zlabel('Z')
    ax2.view_init(23, -151)
    pdf.savefig()
    
    # Représentation des valeurs d'équilibre en fonction de K
    fig, graph3 = plt.subplots()
    graph3.plot(KK,Xeq,KK,Yeq,KK,Xeq_ana,'--',KK,Yeq_ana,'--')
    graph3.set_xlabel('Psupply')
    graph3.set_ylabel('Densités de populations')
    pdf.savefig()
    
    # Représentation de la stabilité en fonction de K
    fig, graph4 = plt.subplots()
    graph4.plot(KK,VP,[gmax1, gmax2],[0,0],'k--')
    graph4.set_xlabel('Psupply')
    graph4.set_ylabel('VPmax')
    pdf.savefig()