"""
Test on parameter effects on phytoplancton concentration
"""

import matplotlib.pyplot as plt
import numpy as npy
import numpy as np
import scienceplots

plt.style.use(['science','no-latex'])
plt.close('all')

from growth_model_2P1Z_v9 import growth_model_2P1Z_v9 

""" Set up """

# Timestep (same than timestep of the growth function : by default is 0.2)
dt = 0.2
# Choose the time at which the simulation end (days)
end_time = 600
# create a time vector from 0 to end_time, seperated by a timestep
time = npy.arange(0,end_time,dt)

# Externe input of nutrients 
Psupply_moy = 1 #0.01
#Psupply constant
Psupply =[Psupply_moy] * len(time)

# !!!! use Psupply_arr instead of Psupply when you call the function to test a range of Psupply values (Psupply_moy must be equal to 1)
Psupply_arr = np.array(Psupply)

""" MODIFY THIS TO CHOOSE THE PARAMETERS AND VALUE RANGE YOU WANT TO TEST """

n = 100
n2 = 9 #to reduce the number of line in the line plots 

#name of the parameters in the model you want to study 
param1 = 'gmax2'
param2 = 'psupp'

# name for plot and pdf
name_param1='GMAX2'
name_param2='PSUPPLY'
name_pdf='outputs/model_V9/test1_gmax2_psupp_v9.pdf'

min_param1 = 0.1
max_param1 = 1.5
l_param1 = np.linspace(min_param1,max_param1,n)

min_param2 = 0.01
max_param2 = 1
l_param2 = np.linspace(min_param2,max_param2,n)

P_1=np.zeros((len(l_param1),len(l_param2)))
P_2=np.zeros((len(l_param1),len(l_param2)))
ratio=np.zeros((len(l_param1),len(l_param2)))

""" loop on i for param1 and j for param2 """ 

# !!!! WARNING : Modify the kwargs name to select the good parameters !!!! #

i = 0
for param1 in l_param1:
    j=0
    for param2 in l_param2:
        [P1, P2, Z, PO4,arg]=growth_model_2P1Z_v9(Psupply_arr*param2,time,gmax2=param1) #A modifier
        P_1[i,j] = P1[len(P1)-1]
        P_2[i,j] = P2[len(P2)-1]
        ratio[i,j] = P1[len(P1)-1]/(P1[len(P1)-1]+P2[len(P2)-1])
        j=j+1
    i=i+1
    print(i)

""" Figure : Evolution of biomasses with variation of parameter(s) """

Var=ratio
Var1=P_1
Var2=P_2

from matplotlib.backends.backend_pdf import PdfPages
with PdfPages(name_pdf) as pdf:
    
    # Figure 1 : P function of param1 or param2 only
    
    #P1
    plt.figure(1)
    plt.plot(l_param1[0:n2], Var1[0:n2,0:n2])
    plt.grid()
    plt.xlabel(name_param1)
    plt.ylabel('P1')
    L_param2 = [round(num, 2) for num in l_param2]
    legend = plt.legend(L_param2,frameon=True,fontsize=6,loc='center left')
    legend.set_title(name_param2)
    pdf.savefig()
    
    #P2
    plt.figure(2)
    plt.plot(l_param1[0:n2], Var2[0:n2,0:n2])
    plt.grid()
    plt.xlabel(name_param1)
    plt.ylabel('P2')
    L_param2 = [round(num, 2) for num in l_param2]
    legend = plt.legend(L_param2,frameon=True,fontsize=6,loc='center left')
    legend.set_title(name_param2)  
    pdf.savefig()

    
    #P1/(P1+P2)
    plt.figure(3)
    plt.plot(l_param1[0:n2], Var[0:n2,0:n2])
    plt.grid()
    plt.xlabel(name_param1)
    plt.ylabel('P1/(P1+P2)')
    L_param2 = [round(num, 2) for num in l_param2]
    legend = plt.legend(L_param2,frameon=True,fontsize=6,loc='center left')
    legend.set_title(name_param2)
    #plt.ylim([0,1])
    pdf.savefig()
    
    # Figure 2 and 3 : P concentration function of param2 vs param1 
    
    #P1
    plt.figure(4)
    [X,Y] = np.meshgrid(l_param1,l_param2, indexing='ij')
    plt.pcolor(X,Y,Var1)
    plt.xlabel(name_param1)
    plt.ylabel(name_param2)
    plt.colorbar()
    plt.title('P1')
    plt.clim(0, 1)
    pdf.savefig()
    
    #P2
    plt.figure(5)
    [X,Y] = np.meshgrid(l_param1,l_param2, indexing='ij');
    plt.pcolor(X,Y,Var2)
    plt.xlabel(name_param1)
    plt.ylabel(name_param2)
    plt.colorbar()
    plt.title('P2')
    plt.clim(0, 1)
    pdf.savefig()
    
    #P1/(P1+P2)
    plt.figure(6)
    [X,Y] = np.meshgrid(l_param1,l_param2, indexing='ij');
    plt.pcolor(X,Y,Var)
    plt.xlabel(name_param1)
    plt.ylabel(name_param2)
    plt.colorbar()
    plt.title('P1/(P1+P2)')
    plt.clim(0, 1)
    pdf.savefig()