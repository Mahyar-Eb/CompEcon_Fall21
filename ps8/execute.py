#!/usr/bin/env python
# coding: utf-8

# In[29]:


#1.packages
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
from scipy.stats import norm
import scipy.integrate as integrate
#2.Parameters
beta= 0.7
rho = 0.4
mu = 1/(1+0.02)
sigma = 0.6
nu_shocks = 5
eps = np.random.normal(0.5, sigma, size=(nu_shocks))
lowerB_i = 0.0 #Lower bound for Investment
UpperB_i = 8.0 #Upper bound for Investment
NG_i = 500  # Number of grid point
VFtol = 1e-15 
VFdist = 5.0 
VFmaxiter = 500
import functions as fn
#3.computing R_t
R = np.empty(nu_shocks)
R=fn.R_t(nu_shocks,rho,mu,eps)
#4.creating grid point for shocks
R_grid, pi= fn.shock_g(nu_shocks, mu, rho, sigma)
# 5.Creating grid for state space
I_g = fn.ss_g(lowerB_i,UpperB_i, NG_i)
#6.utility function
C = np.zeros((nu_shocks, NG_i, NG_i))
U= fn.utility(C,nu_shocks,NG_i, R_grid,I_g )
#7.Value function
VF=fn.val_fun(U,pi,nu_shocks,NG_i,VFmaxiter,VFdist,VFtol,beta)
# Visualization, value function 
plt.figure()
fig, ax = plt.subplots()
plt.scatter(I_g[1:], VF[0,1:], label='$R$ = ' + str(R_grid[0]))
plt.scatter(I_g[1:], VF[1,1:], label='$R$ = ' + str(R_grid[1]))
plt.scatter(I_g[1:], VF[2,1:], label='$R$ = ' + str(R_grid[2]))
plt.scatter(I_g[1:], VF[3,1:], label='$R$ = ' + str(R_grid[3]))
plt.scatter(I_g[1:], VF[4,1:], label='$R$ = ' + str(R_grid[4]))
legend = ax.legend(loc='lower left', shadow=False)
plt.xlabel('Investment')
plt.ylabel('Value Function')
plt.title('Value Function as a function of Investment')
plt.show()
# 9.Policy function visualization
opt_i_0 = I_g[PF[0]] # tomorrow's optimal Investment
opt_i_1 = I_g[PF[1]] 
opt_i_2 = I_g[PF[2]] 
opt_i_3 = I_g[PF[3]] 
opt_i_4 = I_g[PF[4]] 

optC_0 =I_g - opt_i_0 # optimal consumption 
optC_1 =I_g - opt_i_1
optC_2 =I_g - opt_i_2
optC_3 =I_g - opt_i_3
optC_4 =I_g - opt_i_4

plt.figure()
fig, ax = plt.subplots()
ax.plot(I_g[1:], opt_i_0[1:], label='Investment'+'$R$ = ' + str(R_grid[0]))
ax.plot(I_g[1:], opt_i_1[1:], label='Investment'+'$R$ = ' + str(R_grid[1]))
ax.plot(I_g[1:], opt_i_2[1:], label='Investment''$R$ = ' + str(R_grid[2]))
ax.plot(I_g[1:], opt_i_3[1:], label='Investment'+'$R$ = ' + str(R_grid[3]))
ax.plot(I_g[1:], opt_i_4[1:], label='Investment'+'$R$ = ' + str(R_grid[4]))

ax.plot(I_g[1:], I_g[1:], '--', label='45 degree line')
legend = ax.legend(loc='lower right', shadow=False)
plt.xlabel('Investment')
plt.ylabel('Consumption')
plt.title('Policy Function - For Different Levels of Shock')
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:




