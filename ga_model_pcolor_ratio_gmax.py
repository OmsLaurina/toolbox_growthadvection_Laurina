#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 11:39:50 2022

@author: loms
"""
import numpy as np
import scipy
from ypstruct import struct


def ga_model_2P2Z_v1(Nsupply,C_nut,dt,nb,time):
    
    ## -------------- Default parameter
    
    umax_small = 1.9872
    umax_big = 2.7648
    gmax_small = 1.4256
    gmax_big  = 1.4256
    kNreg=13.3										
    kNnew = 13.3												
    kG_small=5											
    kG_big=5											
    mP=0.05														
    mZ=0.005											
    eZ=0.1												
    epsilon=0.75											
    P_small_ini=0.6
    P_big_ini=0.1
    Z_small_ini=0.3
    Z_big_ini=0.2
    
    ## -------------- Initial conditions
    
    Nnew=time*np.nan
    Nreg=time*np.nan 				
    P_small=time*np.nan 			
    P_big=time*np.nan 			
    Z_small=time*np.nan 			
    Z_big=time*np.nan				
    u_big=time*np.nan 			
    u_small=time*np.nan 			
    g_big=time*np.nan 			
    g_small=time*np.nan			
    PP_big=time*np.nan 			
    PP_small=time*np.nan 			
    G_big1=time*np.nan 			
    G_big2=time*np.nan 			
    G_small=time*np.nan 			
    excretion_Zbig=time*np.nan 	
    excretion_Zsmall=time*np.nan 	
    death_Zbig=time*np.nan 		
    death_Pbig=time*np.nan 		
    regeneration_Nreg=time*np.nan
    
    P_small[1]=P_small_ini 
    P_big[1]=P_big_ini 
    Z_small[1]=Z_small_ini 
    Z_big[1]=Z_big_ini 	
    Nnew[1]=C_nut
    Nreg[1]=C_nut 
    
    
    ## -------------- Time
    nb_time=len(time)
    
    ## -------------- Loop on time
    
    for t in range(2,nb_time):
        # growth and grazing rates
        u_big[t-1]=Nnew[t-1]/(kNnew+Nnew[t-1])*umax_big
        u_small[t-1]=Nreg[t-1]/(kNreg+Nreg[t-1])*umax_small
        g_big1=P_big[t-1]/(kG_big+Z_small[t-1]+P_big[t-1])*gmax_big
        g_big2=Z_small[t-1]/(kG_big+Z_small[t-1]+P_big[t-1])*gmax_big
        g_big[t-1]=g_big1+g_big2
        g_small[t-1]=P_small[t-1]/(kG_small+P_small[t-1])*gmax_small
        
        # fluxes
        PP_big[t]=u_big[t-1]*P_big[t-1] 
        PP_small[t]=u_small[t-1]*P_small[t-1] 
        G_big1[t]=g_big1*Z_big[t-1] 
        G_big2[t]=g_big2*Z_big[t-1] 
        G_small[t]=g_small[t-1]*Z_small[t-1]
        death_Zbig[t]=mZ*Z_big[t-1]**2
        death_Pbig[t]=mP*P_big[t-1]
        excretion_Zbig[t]=eZ*Z_big[t-1]
        excretion_Zsmall[t]=eZ*Z_small[t-1]
        regeneration_Nreg[t]=epsilon*excretion_Zbig[t]+excretion_Zsmall[t]	
        maxPP_big=(Nnew[t-1]+Nsupply[t]*dt)/dt 
        maxPP_small=(Nreg[t-1]+regeneration_Nreg[t]*dt)/dt 
        if PP_big[t]>maxPP_big:
            PP_big[t]=maxPP_big
        if PP_small[t]>maxPP_small:
            PP_small[t]=maxPP_small
        
        # carbon-based nutrients and biomass
        Nnew[t]=Nnew[t-1]+Nsupply[t]*dt-PP_big[t]*dt 
        Nreg[t]=Nreg[t-1]+regeneration_Nreg[t]*dt-PP_small[t]*dt 
        P_big[t]=P_big[t-1]+PP_big[t]*dt-G_big1[t]*dt-death_Pbig[t]*dt
        #P_big(P_big<=0)=0
        P_small[t]=P_small[t-1]+PP_small[t]*dt-G_small[t]*dt 
        #P_small(P_small<=0)=0
        Z_small[t]=Z_small[t-1]+G_small[t]*dt-G_big2[t]*dt-excretion_Zsmall[t]*dt 
        #Z_small(Z_small<=0)=0
        Z_big[t]=Z_big[t-1]+G_big1[t]*dt+G_big2[t]*dt-excretion_Zbig[t]*dt-death_Zbig[t]*dt 
        #Z_big(Z_big<=0)=0
        
        
        ## -------------- Ouputs
        
        
        # units=struct('time','days','Nsupply','mmolC m^{-3} d^{-1}',
        #               'P_small','mmolC m^{-3}','P_big','mmolC m^{-3}','Z_small','mmolC m^{-3}','Z_big','mmolC m^{-3}',
        #               'Nnew','mmolC m^{-3}','Nreg','mmolC m^{-3}','Chl','mg m^{-3}','PP','gC m^{-3}/yr',
        #               'u_small','d^{-1}','u_big','d^{-1}','g_small','d^{-1}','g_big','d^{-1}')
        
        # output=np.struct('units',units,'time',time,'Nsupply',Nsupply,
        #               'P_small',P_small,'P_big',P_big,'Z_small',Z_small,'Z_big',Z_big,'Nnew',Nnew,'Nreg',Nreg,
        #               'PP',(PP_big+PP_small)*1E-3*12*365.25,
        #               'u_small',u_small,'u_big',u_big,'g_small',g_small,'g_big',g_big)
        
        return 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
        
        return
    