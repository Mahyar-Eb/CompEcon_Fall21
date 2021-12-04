#!/usr/bin/env python
# coding: utf-8

# In[1]:


#1.Importing packages
import numpy as np 
import ar1_approx as ar
#2.computing R_t
def R_t(nu_shocks,rho,mu,eps ):
    R = np.empty(nu_shocks)
    R[0] = 0.0 + eps[0]
    for t in range(1, nu_shocks):
        R[t] = rho * R[t - 1] + (1 - rho) * mu + eps[t]
    return R[t]
# 4.creating grid point for shocks
def shock_g(nu_shocks, mu, rho, sigma):
    ln_R_grid, pi_t =ar.addacooper(nu_shocks, mu, rho, sigma)
    R_grid = np.exp(ln_R_grid)
    pi = np.transpose(pi_t)
    return R_grid, pi
#5.Creating grid for state space
def ss_g(lowerB_i,UpperB_i, NG_i):
    I_g = np.linspace(lowerB_i,UpperB_i, NG_i)
    return I_g
def utility(C,nu_shocks,NG_i, R_grid,I_g ):
     
    for e in range(nu_shocks):
        for i in range(NG_i):
            for k in range(NG_i):
                C[e,i,k] = R_grid[e] * I_g[i]-I_g[k]

    C[C<=0] = 1e-15 # non-negative constraint for consumption
    
    if (C == 0).any():
        U = 0
    else:
        U = np.log(C)
    U[C<0] = -9999999
    return U

# 7.Value Function Iteration 
'''
------------------------------------------------------------------------
VFtol     = scalar, tolerance required for value function to converge
VFdist    = scalar, distance between last two value functions
VFmaxiter = integer, maximum number of iterations for value function
V         = matrix, the value functions at each iteration
Vmat      = matrix, the value for each possible combination of i and i'and R
Vstore    = matrix, stores V at each iteration 
VFiter    = integer, current iteration number
TV        = vector, the value function after applying the Bellman operator
PF        = vector, indicies of choices of i' for all i and R
VF        = vector, the "true" value function
------------------------------------------------------------------------
'''
def val_fun(U,pi,nu_shocks,NG_i,VFmaxiter,VFdist,VFtol,beta):
     
    V = np.zeros((nu_shocks, NG_i) ) 
    Vmat = np.zeros((nu_shocks,NG_i, NG_i )) 
    Vstore = np.zeros((nu_shocks, NG_i, VFmaxiter)) 
    VFiter = 1 
    while VFdist > VFtol and VFiter < VFmaxiter:
        for i in range(nu_shocks): 
            for j in range(NG_i): 
                    EV = 0
                    for ii in range(nu_shocks): 
                        EV += pi[i, ii] * V[ii, max(j - 1, 0)]   
                    Vmat[i, j] = U[i, j] + beta * EV   
        Vstore[:, VFiter] = V.reshape((nu_shocks,NG_i)) 
        TV = Vmat.max(1) 
        PF = np.argmax(Vmat, axis=1)
        VFdist = (np.absolute(V - TV)).max()  
        V = TV
        VFiter += 1 
        
    if VFiter < VFmaxiter:
        print('Value function converged after this many iterations:', VFiter)
    else:
        print('Value function did not converge')            
    # solution to the functional equation
    VF = V 


# In[ ]:




