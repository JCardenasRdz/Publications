import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  
    
    
def Cp10(t):
    """
    Input:  time t in minutes
    THIS FUNCTION CALCULATES AN AIF WITH A SIMULATED INJECTION TIME OF 10
    SECONDS
    Injection of 10 seconds
    """
    A= 30.0 # mM/min
    B= 1.0  
    C= 4.0  # min^-1
    D= 0.65 # mM
    E= 5.0  # min
    F= 0.04 # min-1
    
    Cp_out = A*(t**B)*np.exp(-t**C)+ D*(1-np.exp(-t*E))* np.exp(-t*F);
    
    return Cp_out
        

def Tofts(pars,time,Cp):
        ktrans=pars[0]
        kep=pars[1]
        
        n_points = len(time)
        c_toi= np.zeros( (n_points,1) )
        
        #c_toi=np.zeros(n_points,1)
        
        return c_toi

    
ctoi = []

for k in np.arange(1,n_points):
    int_t = time[k]
    #print(int_t)
    
    EXPO = []
    CRPEXP = []
    
    for j in np.arange(0,k):
        dummy_time = time[j]
        
        expo   = np.exp(-((kep*(int_t - dummy_time))))
        crpexp = Cp[j]*expo;
        
        
        EXPO.append(expo)
        CRPEXP.append(crpexp)
        
    t2 = time[0:k]
    
    crpexp_integral = np.trapz(t2,np.array(CRPEXP))
    
    ctoi.append(ktrans*crpexp_integral)
    

    
time = np.linspace(0,30*60,100)   
n_points = len(time)

Cp = Cp10(time) 

ktrans = 0.25
ve     = 0.50
kep    = ktrans / ve



d = pd.read_csv('Test Run.csv')
