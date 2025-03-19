"""
Author: Ben Dang - Nguyen Dang 
Email: ppd34@drexel.edu 


This script contains functions describing some neccesary properties of the atmosphere
This script use Standard Atmosphere Model 


"""




import numpy as np 
# I. Pressure function

g0 = 9.80665 #(m/s^2)
R = 8.31413# ( J/mol*K) universal constant 
M = 0.0286944  # (g/mol) - molar mass of dried gas 
mu_0 = 1.77 * 10 ** (-5) # Dynamic density (Pa-s)


# from  0 <= h <= 11000:  T = a1 +b1h 
a1 = 288.04 # Kelvin 
b1 = -0.00649 # (1/m)
# from 11000 < h < 25000: T = 216.54 (K)
a2 = 216.54 
b2 = 0 
# from 25000 < h < 35000: T = a3 + b3h 
a3 = 141.79
b3 = 0.00299
P0 = 101325 # preesure at sea level 

# The range of altitude 
H0 = 0 
H1 = 11000 
H2 = 25000
H3 = 35000  
# Defining the function of Temperature 
def Temperature (h,H0 = H0, H1 = H1, H2 = H2, H3=H3, a1=a1,b1=b1,a2 = a2,b2 = b2,a3=a3,b3=b3):
    """
    Input: h - altitude in meters 
    Output: Temperature of atmosphere in Kelvin (K)
    """
    T = a1
    if 0 <= h <= 11000: 
        T = a1 +b1*h
    elif 11000 < h <= 25000: 
        T = a2 +b2*h
    elif 25000 < h <= 35000: 
        T = a3 +b3*h
    
    return T
# Defining the function pressure 
def Pressure(h,P0 = P0, a=a1,b=b1,M=M,g0=g0,R=R):
    """
    Input: h - altitude in meters 
    Output: pressure of atmosphere in Pascal (Pa)
    """
    P = P0
    if 0 <= h <= 11000: 
        P =P0 * ((a1+b1*h)/a1) ** ( -(M*g0)/(R*b1))
    elif 11000 < h <= 25000:
        P = Pressure (11000) * np.exp(((-M*g0)/(R*a2))* (h-11000))
    elif 25000 < h <= 35000: 
        P = Pressure (25000) * ((a3+b3*h)/a2) ** ( -(M*g0)/(R*b3))
    
    return P 
# Defining the function of density 
def Density (h): 
    """
    Input: h - altitude in meters 
    Output: Density of atmosphere in kg/m^3 
    """
    rho = ((Pressure(h))*(M))/ ((R)*(Temperature(h)))

    return rho 

# Defining the function for viscosity by using Sutherland's equation 

def Viscosity (h,mu_0 = mu_0,T0=a1): 
    """
    Input: h - altitude in meters 
    Output: Viscosity of atmosphere in Pascal.Second (Pa.s)
    """
    v = mu_0 * ((Temperature(h))/(T0)) ** (1.5) * ((Temperature(h)+110.4)/(T0 + 110.4))
    return v 

