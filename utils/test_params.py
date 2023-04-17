"""
Test on parameter effects on phytoplancton concentration
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

from growth_model_2P1Z_v10 import growth_model_2P1Z_v10 

""" Set up """

# Timestep (same than timestep of the growth function : by default is 0.2)
dt = 0.2
# Choose the time at which the simulation end (days)
end_time = 2000
# create a time vector from 0 to end_time, seperated by a timestep
time = npy.arange(0,end_time,dt)

# Externe input of nutrients 

#!! WARNING DEFINE PSUPPLY = 1 IF YOU WANT TO TEST THIS PARAMETER !!
Psupply_moy = 1 #0.01
#Psupply constant
Psupply =[Psupply_moy] * len(time)

# !!!! use Psupply_arr instead of Psupply when you call the function to test a range of Psupply values (Psupply_moy must be equal to 1)!!!
Psupply_arr = np.array(Psupply)

""" MODIFY THIS TO CHOOSE THE PARAMETERS AND VALUE RANGE YOU WANT TO TEST """

number_param = 1

n =  100#number of value tested
n2 = n #to reduce the number of line in the line plots 

#name of the parameters in the model you want to study 
param1 = 'psupp'
param2 = 'kP1'

# name for plot
name_param1='PSUPPLY'
name_param2='kP1'

min_param1 = 0.01
max_param1 = 0.1
l_param1 = np.linspace(min_param1,max_param1,n)

min_param2 = 0.5
max_param2 = 2
l_param2 = np.linspace(min_param2,max_param2,n)

if number_param == 2:
    
    name_pdf=f"../outputs/model_V10/test1_{param1}_{param2}_v10.pdf"
    
    """ loop on i for param1 and j for param2 """ 
    
    P_1=np.zeros((len(l_param1),len(l_param2)))
    P_2=np.zeros((len(l_param1),len(l_param2)))
    ratio=np.zeros((len(l_param1),len(l_param2)))
    
    # !!!! WARNING : Modify the kwargs name to select the good parameters !!!! #

    i = 0
    for param1 in l_param1:
        j=0
        for param2 in l_param2:
            [P1,P2,Z,PO4,Export,u1,u2,g1,g2,arg]=growth_model_2P1Z_v10(Psupply_arr*param1,time,kP1=param2) #A modifier
            P_1[i,j] = P1[len(P1)-1]
            P_2[i,j] = P2[len(P2)-1]
            ratio[i,j] = P1[len(P1)-1]/(P1[len(P1)-1]+P2[len(P2)-1])
            j=j+1
        i=i+1
        print(i)
    
else:
    name_pdf=f"../outputs/model_V10/test1_{param1}_v10.pdf"
    
    'If you want to study only param1'
    
    P_1=np.zeros((len(l_param1)))
    P_2=np.zeros((len(l_param1)))
    ratio=np.zeros((len(l_param1)))
    
    # !!!! WARNING : Modify the kwargs name to select the good parameters !!!! #     
    
    i=0
    for param1 in l_param1:
        [P1,P2,Z,PO4,Export,u1,u2,g1,g2,arg]=growth_model_2P1Z_v10(Psupply_arr*param1,time) #A modifier
        P_1[i] = P1[len(P1)-1]
        P_2[i] = P2[len(P2)-1]
        ratio[i] = P1[len(P1)-1]/(P1[len(P1)-1]+P2[len(P2)-1])
        i=i+1
        print(i)
    
    """ Figure : Evolution of biomasses with variation of parameter(s) """
    
Var=ratio
Var1=P_1
Var2=P_2

from matplotlib.backends.backend_pdf import PdfPages
with PdfPages(name_pdf) as pdf:
    # Figure 1 : P function of param1 or param2 only
    
    if number_param == 2:   
        
        # #P1
        # plt.figure(1)
        # plt.plot(l_param1[0:n2], Var1[0:n2,0:n2])
        # plt.grid()
        # plt.xlabel(name_param1)
        # plt.ylabel('P1')
        # L_param2 = [round(num, 2) for num in l_param2]
        # legend = plt.legend(L_param2,frameon=True,fontsize=6,loc='center left')
        # legend.set_title(name_param2)
        # pdf.savefig()
        
        # #P2
        # plt.figure(2)
        # plt.plot(l_param1[0:n2], Var2[0:n2,0:n2])
        # plt.grid()
        # plt.xlabel(name_param1)
        # plt.ylabel('P2')
        # L_param2 = [round(num, 2) for num in l_param2]
        # legend = plt.legend(L_param2,frameon=True,fontsize=6,loc='center left')
        # legend.set_title(name_param2)  
        # pdf.savefig()
    
        #P1/(P1+P2)
        plt.figure(3)
        cmap = plt.cm.get_cmap('Blues')
        norm = plt.Normalize(vmin=l_param2.min(), vmax=l_param2.max())
        scalar_map = ScalarMappable(norm=norm, cmap=cmap)
        for i in range(n2):
            color_val = scalar_map.to_rgba(l_param2[i])
            plt.plot(l_param1[0:n2], Var[0:n2,i], color=color_val)
        for i in range(0,n2):
            plt.scatter(l_param1[0:n2], Var[0:n2,i], c=P_2[0:n2,i],cmap='Wistia')
        cb =plt.colorbar()
        cb.set_label('[P1]')
        plt.grid()
        plt.xlabel(name_param1)
        plt.ylabel('P1/(P1+P2)')
        L_param2 = [round(num, 2) for num in l_param2]
        legend = plt.legend(L_param2[0:n2],frameon=True,fontsize=6,loc='center left',bbox_to_anchor=(1.5, 0.5))
        legend.set_title(name_param2)  
        pdf.savefig(bbox_inches='tight', pad_inches=0.5)
        
        #Figure 2 and 3 : P concentration function of param2 vs param1 
        
        # #P1
        # plt.figure(4)
        # [X,Y] = np.meshgrid(l_param1,l_param2, indexing='ij')
        # plt.pcolor(X,Y,Var1)
        # plt.xlabel(name_param1)
        # plt.ylabel(name_param2)
        # plt.colorbar()
        # plt.title('P1')
        # plt.clim(0, 1)
        # pdf.savefig()
        
        # #P2
        # plt.figure(5)
        # [X,Y] = np.meshgrid(l_param1,l_param2, indexing='ij');
        # plt.pcolor(X,Y,Var2)
        # plt.xlabel(name_param1)
        # plt.ylabel(name_param2)
        # plt.colorbar()
        # plt.title('P2')
        # plt.clim(0, 1)
        # pdf.savefig()
        
        # #P1/(P1+P2)
        # plt.figure(6)
        # [X,Y] = np.meshgrid(l_param1,l_param2, indexing='ij');
        # plt.pcolor(X,Y,Var)
        # plt.xlabel(name_param1)
        # plt.ylabel(name_param2)
        # plt.colorbar()
        # plt.title('P1/(P1+P2)')
        # plt.clim(0, 1)
        # pdf.savefig()
    
    else:
        
        # #P1
        # plt.figure(1)
        # plt.plot(l_param1[0:n2], Var1[0:n2])
        # plt.grid()
        # plt.xlabel(name_param1)
        # plt.ylabel('P1')
        # pdf.savefig()
        
        #P2
        # plt.figure(2)
        # plt.plot(l_param1[0:n2], Var2[0:n2])
        # plt.grid()
        # plt.xlabel(name_param1)
        # plt.ylabel('P2')
        # pdf.savefig()
        
        #P1/(P1+P2)
        plt.figure(3)
        plt.plot(l_param1[0:n2], Var[0:n2], color='black')
        plt.scatter(l_param1[0:n2], Var[0:n2], c=P_2,cmap='Wistia')
        cb =plt.colorbar()
        cb.set_label('[P1]')
        plt.grid()
        plt.xlabel(name_param1)
        plt.ylabel('P1/(P1+P2)')
        pdf.savefig()
    