"""
Plot the monod function and holling tyoe II for growth and predation rates
"""
import numpy as np
import matplotlib.pyplot as plt

def f_monod(PO4,u_max,kP,color):
    PO4 = np.array(PO4)
    mu = (PO4/(PO4+kP))*u_max
    plt.scatter(PO4,mu,5,color=color)
    return mu

def f_hollingII(P,sum_P,g_max,kZ,color):
    P = np.array(P)
    sum_P = np.array(sum_P)
    g = (P/(sum_P+kZ))*g_max
    plt.scatter(P,g,5,color=color)
    return g