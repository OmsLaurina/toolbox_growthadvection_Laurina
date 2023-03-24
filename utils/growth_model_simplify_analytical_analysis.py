#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simplify model : analytical analysis

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from scipy import integrate
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
gmax1=0.3     # maximum grazing rates of Z on P1 [d**{-1}]
gmax2=0.6     # maximum grazing rates of Z on P2 [d**{-1}]
mP=0		  # P2 mortality rate (default 0 ie no P2 sinking) [d**{-1}]
mZ=0.05 	  # Z2 quadratic mortality rate [mmolC**{-1} m**{3} d**{-1}]
eZ=1          # zoo excretion rate (Z1 and Z2) [d**{-1}]
gamma1=0.6    # conversion factor from P1 to Z [/]
gamma2=0.6    # conversion factor from P2 to Z [/]
epsilon=0.25  # fraction of Z excretion that is available as regenerated PO4 [/]
Psupply = 0.01

def growth_model_simplify(t,u):

    global umax1, umax2, gmax1, gmax2, eZ, epsilon, mZ, gamma1, gamma2, Psupply
    
    p1=u[0];p2=u[1];z=u[2];po4=u[3]

    
    dp1 = po4*umax1*p1-p1*gmax1*z
    dp2 = po4*umax2*p2-p2*gmax2*z
    dz = gamma1*p1*gmax1*z+gamma2*p2*gmax2*z-eZ*((1-gamma1)*p1*gmax1*z+(1-gamma2)*p2*gmax2*z)-mZ*z**2
    dpo4 = Psupply+epsilon*eZ*((1-gamma1)*p1*gmax1*z+(1-gamma2)*p2*gmax2*z)-(po4*umax1*p1+po4*umax2*p2)
    
    return [dp1, dp2, dz, dpo4]

CI=[p10,p20,z0,po40]
sol=solve_ivp(growth_model_simplify,[0,end_time],CI,method='RK45',t_eval=time)
t=sol.t
p1=sol.y[0,:];p2=sol.y[1,:];z=sol.y[2,:];po4=sol.y[3,:]

""" Equilibre du modèle calculés analytiquement """

# 1/ P1 = 0 
print('Equilibrium values for P1 =0 :')
PO4eq1 = Psupply*(gmax2*(gamma2-eZ*(1-gamma2))*(epsilon*eZ*(1-gamma2)+1))/mZ
PO4eq1 = np.array(PO4eq1)
print('PO4 =', PO4eq1)
Zeq1 = (umax2/gmax2)*PO4eq1
Zeq1 = np.array(Zeq1)
print('Z =', Zeq1)
P2eq = mZ*Zeq1/(gamma2*gmax2-eZ*(1-gamma2*gmax2))
P2eq = np.array(P2eq)
print('P2 =',P2eq)

EQ = [0,P2eq, Zeq1, PO4eq1]

# # Matrice Jacobienne
# # J1
# j11 = PO4eq*umax1-gmax1*Zeq;
# j12 = 0;
# j13 = 0;
# j14 = 0;
# j21 = 0;
# j22 = PO4eq*umax2-gmax2*Zeq;
# j23 = PO4eq*umax2*P2eq-gmax2*P2eq;
# j24 = umax2*P2eq-P2eq*gmax2*Zeq;
# j31 = gamma1*gmax1*Zeq+gamma2*P2eq*gmax2*Zeq-eZ*((1-gamma1)*gmax1*Zeq+(1-gamma2)*P2eq*gmax2*Zeq)-mZ*Zeq**2;
# j32 = gamma2*gmax2*Zeq-eZ*(1-gamma2)*gmax2-mZ*Zeq**2;
# j33 = P2eq*gmax2*gamma2-eZ*(1-gamma2)*gmax2-2*mZ*Zeq;
# j34 = 0;
# j41 = Psupply+epsilon*eZ*((1-gamma1)*gmax1*Zeq+(1-gamma2)*P2eq*gmax2*Zeq)-(PO4eq*umax1+PO4eq*P2eq*umax2);
# j42 = Psupply+epsilon*eZ*((1-gamma2)*gmax2*Zeq)-PO4eq*umax2;
# j43 = Psupply+epsilon*eZ*((1-gamma2)*P2eq*gmax2)-PO4eq*umax2*P2eq;
# j44 = Psupply+epsilon*eZ*((1-gamma2)*P2eq*gmax2*Zeq)-umax2*P2eq;
# J1 = [[j11,j12,j13,j14], [j21,j22,j23,j24], [j31,j32,j33,j34], [j41,j42,j43,j44]]
# J1 = np.array(J1)
# vp=np.linalg.eig(J1)
# det = np.linalg.det(J1)
# print('Le determinat de J est :', det)

# 2/ P2 = 0
print('Equilibrium values for P2 =0 :')
PO4eq2 = Psupply*(gmax1*(gamma1-eZ*(1-gamma1))*(epsilon*eZ*(1-gamma1)+1))/mZ
PO4eq2 = np.array(PO4eq2)
print('PO4 =',PO4eq2)
Zeq2 = (umax1/gmax1)*PO4eq2
Zeq2 = np.array(Zeq2)
print('Z =',Zeq2)
P1eq = mZ*Zeq2/(gamma1*gmax1-eZ*(1-gamma1*gmax1))
P1eq = np.array(P1eq)
print('P1 =',P1eq)

# #J2
# j11 = PO4eq*umax1-gmax1*Zeq;
# j12 = 0;
# j13 = PO4eq*umax1*P1eq-gmax1*P1eq;
# j14 = umax1*P1eq-P1eq*gmax1*Zeq;
# j21 = 0;
# j22 = PO4eq*umax2-gmax2*Zeq;
# j23 = 0
# j24 = 0
# j31 = gamma1*gmax1*Zeq-eZ*(1-gamma1)*gmax1*Zeq-mZ*Zeq**2;
# j32 = gamma2*gmax2*Zeq+gamma1*P1eq*gmax1*Zeq-eZ*((1-gamma2)*gmax2*Zeq+(1-gamma2)*P1eq*gmax1*Zeq)-mZ*Zeq**2;
# j33 = P1eq*gmax1*gamma1-eZ*(1-gamma1)*gmax1*P1eq-2*mZ*Zeq;
# j34 = 0;
# j41 = Psupply+epsilon*eZ*(1-gamma1)*gmax1*Zeq-(PO4eq*umax1);
# j42 = Psupply+epsilon*eZ*((1-gamma2)*gmax2*Zeq+(1-gamma1)*P1eq*gmax1*Zeq)-(PO4eq*umax2+PO4eq*P1eq*umax1);
# j43 = Psupply+epsilon*eZ*((1-gamma1)*P1eq*gmax1)-PO4eq*umax1*P1eq;
# j44 = Psupply+epsilon*eZ*((1-gamma1)*P1eq*gmax1*Zeq)-umax1*P1eq;
# J2 = [[j11,j12,j13,j14], [j21,j22,j23,j24], [j31,j32,j33,j34], [j41,j42,j43,j44]]
# J2 = np.array(J2)
# vp=np.linalg.eig(J1)
# det = np.linalg.det(J1)
# print('Le determinat de J est :', det)


""" Figures """

# Temporal evolution
plt.figure(1)
plt.plot(t, p1, label='P1', color="chartreuse")
plt.plot(t, p2, label='P2', color="green")
plt.plot(t, z, label='Z', color="aqua")
plt.plot(t,po4, label='PO4', color="magenta")
plt.xlabel('Time')
plt.ylabel('mmolC.m^-3')
plt.legend(frameon=True) 
plt.title('model outputs (concentration over time)')


#P1
plt.figure(2)
ax = plt.axes(projection ='3d')
ax.plot(p1,po4,z)
ax.set_xlabel('P1')
ax.set_ylabel('PO4')
ax.set_zlabel('Z')
ax.view_init(23, -151)

#P2
plt.figure(3)
ax2 = plt.axes(projection ='3d')
ax2.plot(p2,po4,z)
ax2.set_xlabel('P2')
ax2.set_ylabel('PO4')
ax2.set_zlabel('Z')
ax2.view_init(23, -151)

    

