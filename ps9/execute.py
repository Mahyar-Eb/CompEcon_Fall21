import numpy as np
import matplotlib.pyplot as plt
import SS
#parameters
S=90
#chi_n = np.ones(S) it was vector of ones so it won't effect the results
alpha = 0.35
delta = 0.125
A = 1
beta = 0.84
sigma = 2.5
l_tilda= 1.0
b= 0.501
nu= 1.554
#initial guesses
r_guess = 0.05
#solve for steady state 
r_ss, w_ss, c, n, b_sp ,success = SS.solver(r_guess, beta,alpha, delta, A,sigma, l_tilda, b, nu, S)
print(r_ss, w_ss, c, n, b_sp)
#plot for saving
plt.figure()
plt.plot(np.arange(90),b_sp,color='red')
plt.xlabel('age')
plt.ylabel('saving')
plt.title('Steady-State distribution of saving')
plt.savefig('saving.png')
#plot for consumption
plt.figure()
plt.plot(np.arange(90),c, color='green')
plt.xlabel('age')
plt.ylabel('consumption')
plt.title('Steady-State distribution of consumption')
plt.savefig('consumption.png')
#plot for labor suply
plt.figure()
plt.plot(np.arange(90),n, color='purple')
plt.xlabel('age')
plt.ylabel('labor supply')
plt.title('Steady-State distribution of Labor Supply')
plt.savefig('labor.png')

