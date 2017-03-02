import numpy as np
from scipy.optimize import root
import numpy as np
from sympy.solvers import solvers
from sympy import Symbol
from sympy import nsolve
from sympy import cos,tanh
from sympy import *
import scipy
from sympy.abc import x,y,m,t
import sympy.utilities.lambdify
from numpy import random
import matplotlib.pyplot as plt


n_pairs = 100
grid_size = 30

# generate a grid with a bias
def grid_generate():
    global m_M_pairs
    m_M_pairs = np.zeros(shape=(2,n_pairs))

    for j in range(0,n_pairs): #n_pairs is the number of time pairs
        N = grid_size**2
        x = 0.8
        m = float(np.random.randint(N/10,size=1)) #arbatry N/n to ensure m/N is small


        s=np.random.choice([0,1],size=(N,),p=[(1-m/N),m/N]) # make grid with s_i= 1 with prob m/N
        #print m,np.sum(s)
        s_1=np.squeeze(np.zeros(shape=(1,N)))


        for i in range(0,N,1):
            if s[i,]==1:
                s_1[i,] = np.random.choice([0,1],p=[(1-x),x])  # s_i detected correctly with prob m_M_pairs[1,i] to give s^1 the underdetected version of s
            else:
                s_1[i,]=0

        M=np.sum(s_1)

        m_M_pairs[0,j]=m
        m_M_pairs[1,j]=float(M)

#differential equations d/dT and d/dm* for ML calculation ======================================================================================


m_zero= 800 # mean of gaussian prior on m*
m_sigma=20 # sigma of gaussian pror on m*
T_zero = 10
T_sigma=2
#m*=(m_zero-(m_sigma**2/t_sigma**2)*m1*(t-T_zero)/t)/(1-(m_sigma**2/t_sigma**2)*(t-T_zero)/t)
M1=679
m1=588


m_bar=m_zero - m_sigma*((M1/(0.5+0.5*tanh(m1-((m_zero-(m_sigma**2/T_sigma**2)*m1*(t-T_zero)/t)/(1-(m_sigma**2/T_sigma**2)*(t-T_zero)/t))/t)))-m1)*(m1 - ((m_zero-(m_sigma**2/T_sigma**2)*m1*(t-T_zero)/t)/(1-(m_sigma**2/T_sigma**2)*(t-T_zero)/t)))*0.5/(t**2)*(1-tanh(m1-((m_zero-(m_sigma**2/T_sigma**2)*m1*(t-T_zero)/t)/(1-(m_sigma**2/T_sigma**2)*(t-T_zero)/t))/t)**2)

print m_bar
m_bar_func=lambdify((t),m_bar,"numpy")
#===================================================================================================


T_old = 10
T_new=2
gamma = 1000 # step size
iterations=0
check=0
printData=true
maxIterations=1000
precision = 0.00001
print m_bar_func(T_old)


while True:
    temp_T= T_old - gamma*m_bar_func(T_old)

    iterations += 1
    if iterations > maxIterations:
        print("Too much iterations. Adjust alpha and make sure that the function be convex!")
        printData = False
        break

    #If the value of theta changes less of a certain amount, our goal is met.
    if abs(temp_T-T_old) < precision:
        break

    #Simultaneous update
    T_old=temp_T
if printData:
    print("The function "+str(expr)+" converges to a minimum")
    print("Number of iterations:",iterations)
    print("theta (x0) =",temp_T)





"""
# iteratively search through m* and T to optimise their values  ====================================================
def optimise(): # search algorithm to find T and m8 iteratively

    dm=0.1 #initial change in m* and T
    dT=0.1
    m_new = 100.0 # initial starting point of m* and T
    T_new = 100.0
    iterations = 1000
    #for steps in range (1,iterations,1):
    #new_func = T_pydif_ML_func(m_new,T_new) # first iteration with initial parameters

    global sim_result
    sim_result = np.zeros(shape = (3,iterations))

    for i in range(1,iterations,1):
            j = np.random.rand()
            if j >= 0.5:

                #m* update
                m_old_func = m_pydif_ML_func(m_new,T_new)
                m_old=m_new

                m_new = m_old + dm
                m_new_func = m_pydif_ML_func(m_new,T_new) # evaluation at the new parameters

                #m_new = m_old
                #change = dm*(new_func-old_func)*(np.random.rand())
                #print change
                #m_new=m_new-change
                #print m_new
                m_pydif_ML_func(m_new,T_new)
                if ((new_func<0)==(old_func<0)): # make sure points are not eitehr side of the root

                    if (abs(new_func)) >= (abs(old_func)): # if the change is a bad one (not closer to the min)
                        m_new = m_old                      # dont accept new m* and change by a smaller step
                        dm=-0.2*dm
                        #print old_func, new_func, m_new, "m_new", T_new, "T_new"
                    else:
                        dm = 1.05 * dm
                        #print old_func, new_func, m_new, "m_new", T_new, "T_new"

                else: # if points lie either side search in between
                    m_new=m_old
                    dm = dm*0.2

                    #print old_func, new_func, m_new, "m_new", T_new, "T_new"

            if j < 0.5:
                #T update
                old_func = T_pydif_ML_func(m_new,T_new)
                T_old=T_new
                T_new = T_old + dT
                new_func = T_pydif_ML_func(m_new,T_new) # evaluation at the new parameters

                #T_new = T_old
                #change = dT*(new_func-old_func)*(np.random.rand())
                #print change
                #t_new=T_new-change
                #print T_new
                T_pydif_ML_func(m_new,T_new)
                if ((new_func<0)==(old_func<0)): # make sure points are not eitehr side of the root

                    if (abs(new_func)) >= (abs(old_func)): # if the change is a bad one (not closer to the min)
                        T_new = T_old                      # dont accept new m* and change by a smaller step
                        dT=-0.2*dT
                        #print old_func, new_func, m_new, "m_new", T_new, "T_new"
                    else:
                        dT = 1.05 * dT
                        #print old_func, new_func, m_new, "m_new", T_new, "T_new"

                else: # if points lie either side search in between
                    T_new=T_old
                    dT = dT*0.2

                    #print old_func, new_func, m_new, "m_new", T_new, "T_new"

            sim_result[0,i-1]=new_func
            sim_result[1,i-1]=m_new
            sim_result[2,i-1]=T_new

plt.plot(sim_result[0,:],sim_result[1,:])
plt.show()
"""
