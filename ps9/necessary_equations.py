#1.import packages
import scipy.optimize as opt
import numpy as np
#Household functions 
def get_c(first_c,r, beta, sigma, p ):
    '''
    finding consumption given initial guess of first_c

    '''
    c = np.zeros(p)
    c[0] = first_c
    c_sp = first_c
    for t in np.arange(p):
        if t<p-1:
                c[t + 1] = c_sp * (beta * (1 + r[t + 1])) ** (1 / sigma)
                c_sp = c[t + 1]
                t += 1
        else:
            break
    return c

def mu_c_fun(c,sigma):
    '''
    Marginal utility of consumption

    '''
    eps = 1e-5
    muc = c ** -sigma
    mu1 = (-sigma) * eps ** (-sigma - 1)
    mu2 = eps ** (-sigma) - mu1 * eps
    c_constr= c < eps
    muc[c_constr] = mu1 * c[c_constr] + mu2
    return muc
def mu_n_fun(n,b, l_tilda, nu):
    '''
    find marginal utility for laber supply, eq.4.9 LHS
    chi = vector, values for χ_s, ones for each period
    b=current period wealth
    l_tilda = scalar, per-period time endowment
    nu= scalar
    I found out that the problem for finding optimal root was the labor disutility
    function I had to define the conditions to satisfy the upper bound and lower bond
    conditions. for this purpose I used the method from the refrence below. I redefine the code 
    as you mentioned but still root finder didn't run so I thought may be it's due to the 
    disutility function.thats why I used my last codes again.
    refrence:
    https://github.com/OpenSourceEcon/LaborCalibrate
    
    '''
    eps_lowerb = 1e-6
    eps_upperb = l_tilda - eps_lowerb 
    mun = (b/l_tilda)*((np.abs(n)/l_tilda)**(nu-1))*((1-((np.abs(n)/l_tilda)** nu))**((1-nu)/nu))
    mu1_l= (b * (l_tilda ** (-nu)) * (nu - 1) * (eps_lowerb  ** (nu - 2)) * \
         ((1 - ((eps_lowerb  / l_tilda) ** nu)) ** ((1 - nu) / nu)) * \
         (1 + ((eps_lowerb  / l_tilda) ** nu) * ((1 - ((eps_lowerb  / l_tilda) ** nu)) ** (-1))))
    mu2_l= ((b / l_tilda) * ((eps_lowerb / l_tilda) ** (nu - 1)) * \
         ((1 - ((eps_lowerb / l_tilda) ** nu)) ** ((1 - nu) / nu)) - (mu1_l * eps_lowerb))
    mu1_up = (b * (l_tilda ** (-nu)) * (nu - 1) * (eps_upperb ** (nu - 2)) * \
         ((1 - ((eps_upperb / l_tilda) ** nu)) ** ((1 - nu) / nu)) * \
         (1 + ((eps_upperb / l_tilda) ** nu) * ((1 - ((eps_upperb / l_tilda) ** nu)) ** (-1))))
    mu2_up = ((b / l_tilda) * ((eps_upperb / l_tilda) ** (nu - 1)) * \
         ((1 - ((eps_upperb / l_tilda) ** nu)) ** ((1 - nu) / nu)) - (mu1_up * eps_upperb))
    n_lower_constr = n < eps_lowerb
    n_upper_constr = n > eps_upperb

    mun[n_lower_constr] = mu1_l * n[n_lower_constr] + mu2_l
    mun[n_upper_constr] = mu1_up * n[n_upper_constr] + mu2_up
    return mun
def get_bsp(c, n ,r, w, p, b_s=0.0):
    
    '''
    calculating lifetime savings, given consumption and labor decisions
    using HH bc, equation

    '''
    b_sp1 = np.zeros(p)
    b_sp1[0]=b_s
    for t in range(p):
        if t<p-1:
            b_sp1[t + 1] = (1 + r[t]) * b_s + w[t] * n[t] - c[t]
            b_s = b_sp1[t + 1]
            t += 1 
        else:
            break
       
    return b_sp1
def get_ee(n,c,w, sigma, l_tilda, b, nu):
    
    '''
    function for solving n given r,w,form households function,equation4.9

    '''
    muc= mu_c_fun(c,sigma)
    mun= mu_n_fun(n,b, l_tilda, nu)
    euler_error = w * muc - mun
    return euler_error
def get_n(c,sigma, l_tilda, b, nu, w):
    '''
    # using initial n guess and euler equation to calculate for n

    ''' 
    n_guess= 0.4 * np.ones(90)
    sol_n = opt.root(get_ee,n_guess, args=(w,c,sigma, l_tilda, b, nu), method = 'lm')
    if sol_n.success:
        n = sol_n.x
    else:
        raise ValueError("Can't find optimal root(labor supply satisfied the euler equation)")
    return n
def get_bl(first_c, r, w, beta, sigma, l_tilda, b, nu, p ): 
    
    ''' 
    last period saving, using hh's BC only for last period to find the 
    the optiomal consumption condition on b_l=0
    b_l= last period saving'''
    # function for last-period savings, given intial guess c1
    c = get_c(first_c, r, beta, sigma, p)
    n = get_n(c,sigma, l_tilda, b, nu, w)
    b_sp1 = get_bsp(c,n ,r, w , p,b_s=0)
    b_l = (1 + r[-1]) * b_sp1[-1] + w[-1] * n[-1] - c[-1]
    return b_l  
# Firms functions
def get_r(K,L,alpha, delta, A):
    '''
    --------------------------------------------------------------------
    Compute the interest rate from the firm's FOC, equation 4.13.chapter 4
    K= scalar, aggregate capital supplied from part 3
    L= scalar, aggregate labor supplied from part 3
    A = scalar, total factor priductivity
    alpha = scalar, captal share of income, α ∈ (0,1)
    delta = scalar, depreciation rate of capital, δ ∈ (0,1)
    r = scalar, steady state interest rate
    '''
    r = alpha *A *(L/K)**(1-alpha)-delta
    return r

def get_w(r,alpha, delta, A):
    '''
    --------------------------------------------------------------------
    solve for wage from the firm's FOC, equation 4.14.chapter 4
    K= scalar, aggregate capital supplied from part 3
    L= scalar, aggregate labor supplied from part 3
    A = scalar, total factor priductivity
    alpha = scalar, captal share of income, α ∈ (0,1)
    delta = scalar, depreciation rate of capital, δ ∈ (0,1)
    w = scalar, steady state wage
    ''' 
    w = (1 - alpha) * A * ((alpha * A) / (r + delta)) ** (alpha / (1 - alpha))
    return w

# function to compute aggregate labor supplied
def get_L(n):
    '''
    --------------------------------------------------------------------
    Computing aggregate labor supplied by using market clearing, equation
    4.15 chapter 4
    n = vector of labor supplied each period
    ''' 
    L = n.sum()
    return L
#3. function to compute aggregate capital supplied
def get_K(b_sp):
    '''
    --------------------------------------------------------------------
    Computing aggregate capital supplied by using market clearing, equation
    4.16 chapter 4
    b = vector of capital supplied each period
    '''
    K = b_sp.sum()
    return K








