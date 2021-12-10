import numpy as np
import scipy.optimize as opt
import necessary_equations as ne
#SS algorithm
    # 1.guess r -> w
    # 2.using r and w from 1, calculate c which satisfies the b_l = 0 -> first_c 
    # 3.using first_c solve for c using get_c function
    # 4.using first_c solve for n using get_n function 
    # 5.solve for b_sp1 using c and n from 3 and 4 
    # 6.solve for c,n,and b to get K and L
    # 7.using K, L in get_r -> r'
    # 8.check if r' = guess of r
    # 6.loop again if not

def solver(r_guess, beta,alpha, delta, A,sigma, l_tilda, b, nu, S):
    # Find Steady State
    max_iter = 500
    iter = 0
    tol = 1e-8
    xi = 0.7
    dist = 8
    #initial guess for r
    while (dist > tol) and (iter < max_iter):
        r = r_guess * np.ones(S)
        w = ne.get_w(r, alpha, delta, A) * np.ones(S)
    #  solve for c which make last period saving equal to zero   
        c_guess= 0.2 
        sol_optc = opt.root(ne.get_bl,c_guess, args=(r, w, beta ,sigma ,l_tilda , b, nu , S))
        if sol_optc.success :
            first_c = sol_optc.x
        else:
            raise ValueError("cannot find optimum root, it doesn't satisfy b_l = 0")
        
        c = ne.get_c(first_c, r, beta, sigma, S)
        n = ne.get_n(c, sigma, l_tilda, b, nu, w)
        b_bsp = ne.get_bsp(c, n, r, w, S)
        K = ne.get_K(b_bsp)
        L = ne.get_L(n)
        r_prim = ne.get_r(K, L, A, alpha, delta)
        dist = ((r_prim - r) ** 2).max()
        iter += 1
        r = xi * r + (1 - xi) * r_prim
    success = iter<max_iter
    return r, w, c, n, b_bsp ,success
