# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""
import matplotlib.pyplot as plt
import numpy as np

def monod(N,u_max,kN):
    
    N = np.array(N)
    
    u = N/(N+kN)
    
   
    if np.size(u) >1:
        plt.scatter(N,u)
        
    return u
