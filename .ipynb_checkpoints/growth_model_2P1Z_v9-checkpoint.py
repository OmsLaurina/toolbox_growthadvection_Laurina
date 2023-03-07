"""
        ----------- Growth_MODEL2P1Z_FROMPSUPPLY [mmolC/m^3] -----------

 Growth-advection model by Messié & Chavez (2017), adapted to the Mediterranean sea.

 Call it as : [P1,P2,Z,PO4]=growth_model_2P1Z_v9(Psupply, time,**kwargs)
 
 Required inputs:
 	'Psupply' expressed in mmolC/m3/d, all units are carbon-based.
    'time' as npy.arange(0,end_time,dt) where dt is the time step and end_time the number of day

 Optional inputs:
     **kwargs parameters of models, initial conditions, timestep. Set by default here, but you can modify it as :
      [P1,P2,Z,PO4]=growth_model_2P1Z_v9(Psupply, time, umax1=0.2) for exemple

 Monique Messié, 2021 for public version
 Reference: Messié, M., & Chavez, F. P. (2017). PO4rient supply, surface currents, and plankton dynamics predict zooplankton hotspots 
					in coastal upwelling systems. Geophysical Research Letters, 44(17), 8979-8986, https://doi.org/10.1002/2017GL074322
                    
 Differences with Messié and Chavez (2021):
	   Nutrient supply is Phosphate instead Nitrate
	   One nutrient compartiment
       One zooplankton compartiment
       Assimilation coefficient gamma = %  of conversion from P to Z
       
  P1 = Small phytoplankton (SYNECO, PICO) < 2um
  P2 = Big phytoplankton (MICRO) >20um
  
"""

def growth_model_2P1Z_v9(Psupply,time,**kwargs):
    
    global dt, umax1, umax2, gmax1, gmax2, kP1, kP2, kZ, eZ, epsilon, mZ, gamma1, gamma2, u1, u2, g1, g2, PP1, PP2, G1, G2, P1, P2, Z, PO4, P01, P02, Z0, PO40

    """" definition of model variables """
    
    ##Creating empty vectors
    
    # growth rates [d^-1]
    u1=[]
    u2=[]		
    # grazing rates [d^-1]
    g1=[]
    g2=[]				
    #Primary production [mmolC^{-1} m^{3} d^{-1}]
    PP1=[]
    PP2=[]
    #Grazing [mmolC^{-1} m^{3} d^{-1}]
    G1=[]
    G2=[]
    # biomass of Phytoplankton and Zooplankton [mmolC/m³]
    P1=[] 			
    P2=[] 						
    Z = []
    #Phosphate concentration [mmolC/m³]
    PO4 = []
    #exc and regeneration [mmolC^{-1} m^{3} d^{-1}]
    exc=[]
    reg=[]
    p_fec=[]
    #death [mmolC^{-1} m^{3} d^{-1}]
    d_Z=[]
    #export [mmolC^{-1} m^{3} d^{-1}]
    Export=[]
    
    """ default value of model parameters and initial values """
    
    defaultKwargs = {
        'umax1':1.9872, # maximum growth rates of P1 [d^{-1}]
        'umax2':2.7648, # maximum growth rates of P2 [d^{-1}]
        'gmax1':0.5,    # maximum grazing rates of Z on P1 [d^{-1}]
        'gmax2':0.4,    # maximum grazing rates of Z on P2 [d^{-1}]
        'kP1':1,        # half-saturation constant for P1 on PO4 [mmolC m^{-3}]
        'kP2':2,	    # half-saturation constant for P2 on PO4 [mmolC m^{-3}]
        'kZ':5,			# half-saturation constant for Z on P [mmolC m^{-3}]
        'mP':0,			# P2 mortality rate (default 0 ie no P2 sinking) [d^{-1}]
        'mZ':0.07,		# Z2 quadratic mortality rate [mmolC^{-1} m^{3} d^{-1}]
        'eZ':0.1,	    # zoo excretion rate (Z1 and Z2) [d^{-1}]
        'gamma1':0.7,   # conversion factor from P1 to Z [/]
        'gamma2':0.7,   # conversion factor from P2 to Z [/]
        'epsilon':0.25, # fraction of Z excretion that is available as regenerated PO4 [/]
        'P1_ini':0.6,   # initial biomass of P1 [mmolC m^{-3}]
        'P2_ini':0.1,   # initial biomass of P2 [mmolC m^{-3}]
        'Z_ini':0.6,    # initial biomass of Z [mmolC m^{-3}]
        'PO4_ini':0.5,  # initial biomass of PO4 [mmolC m^{-3}]
        'DT':0.2        # timestep of discretisation  [day]
        }
    
    arg = { **defaultKwargs, **kwargs}
    
    dt=arg['DT']
    
    # Initial conditions at time=0			
    P1.append(arg['P1_ini'])
    P2.append(arg['P2_ini'])
    Z.append(arg['Z_ini'])
    PO4.append(arg['PO4_ini'])
    
    
    nb_time=len(time)
    
    for t in range(0,nb_time-1):

        # growth rates (monod function)
        u1next=PO4[t-1]/(arg['kP1']+PO4[t-1])*arg['umax1']
        u2next=PO4[t-1]/(arg['kP2']+PO4[t-1])*arg['umax2']
       
        u1.append(u1next)
        u2.append(u2next)
        
        # functional response = grazing rates (holling type II)
        g1next=P1[t-1]/(arg['kZ']+P1[t-1]+P2[t-1])*arg['gmax1']
        g2next=P2[t-1]/(arg['kZ']+P1[t-1]+P2[t-1])*arg['gmax2']
       
        g1.append(g1next)
        g2.append(g2next)
   
   	    # fluxes 
        PP1next=u1[t]*P1[t-1]
        PP2next=u2[t]*P2[t-1]
       
        PP1.append(PP1next)
        PP2.append(PP2next)
       
        G1_next=arg['gamma1']*g1[t]*Z[t-1]
        G2_next=arg['gamma2']*g2[t]*Z[t-1]
        d_Znext=arg['mZ']*Z[t-1]**2
        excnext=arg['eZ']*((1-arg['gamma1'])*g1[t]*Z[t-1]+(1-arg['gamma2'])*g2[t]*Z[t-1]);
       
        G1.append(G1_next)
        G2.append(G2_next)
        d_Z.append(d_Znext)
        exc.append(excnext)
       
        regnext=arg['epsilon']*exc[t]
        p_fecnext=(1-arg['epsilon'])*exc[t]
       
        reg.append(regnext)
        p_fec.append(p_fecnext)
        
        p_fecnext = (1-arg['epsilon'])*exc[t]
        Exportnext =p_fec[t]+d_Z[t]
       
        p_fec.append(p_fecnext)
        Export.append(Exportnext)
       
        #PP1 puise en premier dans le pool de nutriments
        max_available_PO4 = PO4[t-1]+Psupply[t]*dt+reg[t]*dt
       
        if PP1[t]>max_available_PO4/dt:
            PP1[t]=max_available_PO4/dt
            
        max_available_PO4=max_available_PO4-PP1[t]*dt
        if PP2[t]>max_available_PO4/dt:
            PP2[t]=max_available_PO4/dt
            
        
        # biomasses
        P1_next=P1[t-1]+PP1[t]*dt-G1[t]*dt
        if P1_next<=0:
            P1_next=0
        
        P2_next=P2[t-1]+PP2[t]*dt-G2[t]*dt
        if P2_next<=0:
            P2_next=0
            
        Z_next=Z[t-1]+G1[t]*dt+G2[t]*dt-exc[t]*dt-d_Z[t]*dt
        if Z_next<=0:
            Z_next=0
            
        PO4_next=PO4[t-1]+Psupply[t]*dt+reg[t]*dt-PP1[t]*dt-PP2[t]*dt
        if PO4_next<=0:
            PO4_next=0
        
        P1.append(P1_next)
        P2.append(P2_next)
        Z.append(Z_next)
        PO4.append(PO4_next)
        
    return P1,P2,Z,PO4,arg