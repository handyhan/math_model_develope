""" optimization code test"""
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

# From calculation, it is expected that the local minimum occurs at x=9/4

x_old = 30 # The value does not matter as long as abs(x_new - x_old) > precision
y_old = 10
gamma = 0.01 # step size
iterations=0
check=0
printData=true
maxIterations=1000
precision = 0.00001

expr = ((x-4)**2+(y-4)**2 + 25)
diff_x = diff(expr,x)
diff_y = diff(expr,y)



x_equ_func=lambdify((x,y),diff_x,"numpy")
y_equ_func=lambdify((x,y),diff_y,"numpy")


while True:
    temp_x= x_old - gamma*x_equ_func(x_old,y_old)
    temp_y= y_old - gamma*y_equ_func(x_old,y_old)

    iterations += 1
    if iterations > maxIterations:
        print("Too much iterations. Adjust alpha and make sure that the function be convex!")
        printData = False
        break

    #If the value of theta changes less of a certain amount, our goal is met.
    if abs(temp_x-x_old) < precision and abs(temp_y-y_old) < precision:
        break

    #Simultaneous update
    x_old = temp_x
    y_old = temp_y

if printData:
    print("The function "+str(expr)+" converges to a minimum")
    print("Number of iterations:",iterations)
    print("theta (x0) =",temp_x)
    print("theta1 (y0) =",temp_y)



print x_equ_func,y_equ_func

"""
func = np.zeros(shape=(1))
x_func=np.zeros(shape=(1))
print func.shape
while abs(x_new - x_old) > precision:
    j = np.random.rand()
    if j >= 0.5:

        x_new = -gamma * df(x_old)
        func= np.append(func,df(x_old))
        x_func = np.append(x_func,x_old)
    if j < 0.5:
        print func

plt.plot(func,x_func)
plt.show()

print("The local minimum occurs at %d" % x_new)

#optimise()


#plt.plot(sim_result[0,:],sim_result[1,:])
#plt.show()
"""
