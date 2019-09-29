# =============================================================================
# Chebyshev Approximation of functions
# =============================================================================
from sympy import diff;
from sympy import symbols;
import numpy as np;
import matplotlib.pyplot as plt
from matplotlib.pyplot import ylim
from matplotlib.pyplot import xlim
from matplotlib.pyplot import legend
from matplotlib.pyplot import title
from scipy.interpolate import interp1d
from numpy.polynomial.chebyshev import chebfit
from numpy.polynomial.chebyshev import chebval
from numpy.polynomial.chebyshev import chebroots
from numpy.polynomial.chebyshev import chebval2d
import statsmodels.api as sm
from mpl_toolkits.mplot3d import Axes3D
import math as mt
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from sympy import solve

x = np.linspace(-1, 1, num = 20, endpoint = True) #First we define the x space where we are working on
#Now, let's define the three functions for this exercise
y1=np.exp(1/x)
y2=1/(1+25*x**2) #Runge function
y3=(x+abs(x))/2  #Ramp function


#Interpolation with Chebyshev nodes:
def y1(x):
    return np.exp(1/x)
def y2(x):
    return 1/(1+25*x**2)
def y3(x):
    return (x+abs(x))/2 
x = np.linspace(-1, 1, num = 20, endpoint = True)
roots = np.polynomial.chebyshev.chebroots(x)

y1_eval = y1(roots)
y2_eval = y2(roots)
y3_eval = y3(roots)

c1_3_cheb = np.polyfit(roots, y1_eval, 3)
c2_3_cheb = np.polyfit(roots, y2_eval, 3)
c3_3_cheb = np.polyfit(roots, y3_eval, 3)
c1_5_cheb = np.polyfit(roots, y1_eval, 5)
c2_5_cheb = np.polyfit(roots, y2_eval, 5)
c3_5_cheb = np.polyfit(roots, y3_eval, 5)
c1_10_cheb = np.polyfit(roots, y1_eval, 10)
c2_10_cheb = np.polyfit(roots, y2_eval, 10)
c3_10_cheb = np.polyfit(roots, y3_eval, 10)
a1_3_cheb = np.polyval(c1_3_cheb, roots) 
a2_3_cheb = np.polyval(c2_3_cheb, roots) 
a3_3_cheb = np.polyval(c3_3_cheb, roots) 
a1_5_cheb = np.polyval(c1_5_cheb, roots) 
a2_5_cheb = np.polyval(c2_5_cheb, roots) 
a3_5_cheb = np.polyval(c3_5_cheb, roots) 
a1_10_cheb = np.polyval(c1_10_cheb, roots) 
a2_10_cheb = np.polyval(c2_10_cheb, roots) 
a3_10_cheb = np.polyval(c3_10_cheb, roots) 

#Plots of the Chebyshev nodes approximation first function:
plt.plot(roots, y1_eval,'o', label = 'Original function')
plt.plot(roots, a1_3_cheb,'-', label = 'Interpolation order 3')
plt.plot(roots, a1_5_cheb,'--', label = 'Interpolation order 5')
plt.plot(roots, a1_10_cheb,':', label = 'Interpolation order 10')
plt.legend(loc = 'upper right')
plt.xlim(xmin = -1, xmax = 1)
plt.title('Chebyshev monomial interpolations - Function 1', size=10)
plt.ylabel('f(x)', size=9)
plt.xlabel('x', size=9)
plt.show()

#Plots of the Chebyshev nodes approximation runge function:
plt.plot(roots, y2_eval,'o', label = 'Original function')
plt.plot(roots, a2_3_cheb,'-', label = 'Interpolation order 3')
plt.plot(roots, a2_5_cheb,'--', label = 'Interpolation order 5')
plt.plot(roots, a2_10_cheb,':', label = 'Interpolation order 10')
plt.legend(loc = 'upper right')
plt.xlim(xmin = -1, xmax = 1)
plt.title('Chebyshev monomial interpolations - Runge Function ', size=10)
plt.ylabel('f(x)', size=9)
plt.xlabel('x', size=9)
plt.show()

#Plots of the Chebyshev nodes approximation ramp function:
plt.plot(roots, y3_eval,'o', label = 'Original function')
plt.plot(roots, a3_3_cheb,'-', label = 'Interpolation order 3')
plt.plot(roots, a3_5_cheb,'--', label = 'Interpolation order 5')
plt.plot(roots, a3_10_cheb,':', label = 'Interpolation order 10')
plt.legend(loc = 'upper right')
plt.xlim(xmin = -1, xmax = 1)
plt.title('Chebyshev monomial interpolations - Ramp Function ', size=10)
plt.ylabel('f(x)', size=9)
plt.xlabel('x', size=9)
plt.show()



#We compute the associated errors terms Chebyshev interpolation with monomials:
e1_3_cheb = abs(y1_eval-a1_3_cheb)
e2_3_cheb = abs(y2_eval-a2_3_cheb)
e3_3_cheb = abs(y3_eval-a3_3_cheb)
e1_5_cheb = abs(y1_eval-a1_5_cheb)
e2_5_cheb = abs(y2_eval-a2_5_cheb)
e3_5_cheb = abs(y3_eval-a3_5_cheb)
e1_10_cheb = abs(y1_eval-a1_10_cheb)
e2_10_cheb = abs(y2_eval-a2_10_cheb)
e3_10_cheb = abs(y3_eval-a3_10_cheb)

#Now, we do the plot for the error terms for the first function:
plt.plot(roots, e1_3_cheb,'-', label = 'Error of order 3')
plt.plot(roots, e1_5_cheb,'--', label = 'Error of order 5')
plt.plot(roots, e1_10_cheb, ':', label = 'Error of order 10')
plt.legend(loc = 'upper right')
plt.xlim(xmin = -1, xmax = 1)
plt.title('Errors Chebyshev monomial interpolations - Function 1', size=10)
plt.ylabel('f(x)', size=9)
plt.xlabel('x', size=9)
plt.show()


#Now, we do the plot for the error terms for the runge function:
plt.plot(roots, e2_3_cheb,'-', label = 'Error of order 3')
plt.plot(roots, e2_5_cheb,'--', label = 'Error of order 5')
plt.plot(roots, e2_10_cheb, ':', label = 'Error of order 10')
plt.legend(loc = 'upper right')
plt.xlim(xmin = -1, xmax = 1)
plt.title('Errors Chebyshev monomial interpolations - Runge Function', size=10)
plt.ylabel('f(x)', size=9)
plt.xlabel('x', size=9)
plt.show()

#Now, we do the plot for the error terms for the ramp function:
plt.plot(roots, e3_3_cheb,'-', label = 'Error of order 3')
plt.plot(roots, e3_5_cheb,'--', label = 'Error of order 5')
plt.plot(roots, e3_10_cheb, ':', label = 'Error of order 10')
plt.legend(loc = 'upper right')
plt.xlim(xmin = -1, xmax = 1)
plt.title('Errors Chebyshev monomial interpolations - Ramp Function', size=10)
plt.ylabel('f(x)', size=9)
plt.xlabel('x', size=9)
plt.show()

# Chebyshev approximation with Chebyshev polynomial:
#Firstly, we need to create the functions with unknowns
def y1(x):
    return np.exp(1/x)
def y2(x):
    return 1/(1+25*x**2)
def y3(x):
    return (x+abs(x))/2 

vec = np.linspace(-1, 1, num=20, endpoint=True)
roots= np.polynomial.chebyshev.chebroots(vec)

y1 = y1(roots)
y2 = y2(roots)
y3 = y3(roots)
    
# With chebfit we obtain the coefficients of the Chevyshev polynomial, and chebval constructs the polynomial

c13 = np.polynomial.chebyshev.chebfit(roots, y1, 3)
c23 = np.polynomial.chebyshev.chebfit(roots, y2, 3)
c33 = np.polynomial.chebyshev.chebfit(roots, y3, 3)
c15 = np.polynomial.chebyshev.chebfit(roots, y1, 5)
c25 = np.polynomial.chebyshev.chebfit(roots, y2, 5)
c35 = np.polynomial.chebyshev.chebfit(roots, y3, 5)
c110 = np.polynomial.chebyshev.chebfit(roots, y1, 10)
c210 = np.polynomial.chebyshev.chebfit(roots, y2, 10)
c310 = np.polynomial.chebyshev.chebfit(roots, y3, 10)
a13cheb = np.polynomial.chebyshev.chebval(roots, c13)
a23cheb = np.polynomial.chebyshev.chebval(roots, c23)
a33cheb = np.polynomial.chebyshev.chebval(roots, c33)
a15cheb = np.polynomial.chebyshev.chebval(roots, c15)
a25cheb = np.polynomial.chebyshev.chebval(roots, c25)
a35cheb = np.polynomial.chebyshev.chebval(roots, c35)
a110cheb = np.polynomial.chebyshev.chebval(roots, c110)
a210cheb = np.polynomial.chebyshev.chebval(roots, c210)
a310cheb = np.polynomial.chebyshev.chebval(roots, c310)

#Plots of the Chebyshev nodes approximation with Chebyshev polynomial first function:

plt.plot(roots, y1, label = 'Original function')
plt.plot(roots, a13cheb,'-', label = 'Chebyshev order 3')
plt.plot(roots, a15cheb, '--', label = 'Chebyshev order 5')
plt.plot(roots, a110cheb, ':', label = 'Chebyshev order 10')
plt.legend(loc = 'upper right')
plt.xlim(xmin = -1, xmax = 1)
plt.title('Chebyshev polynomial approx. - Function 1', size=10)
plt.ylabel('f(x)', size = 9)
plt.xlabel('x', size = 9)
plt.show()

plt.plot(roots, y2, label = 'Original function')
plt.plot(roots, a23cheb,'-', label = 'Chebyshev order 3')
plt.plot(roots, a25cheb, '--', label = 'Chebyshev order 5')
plt.plot(roots, a210cheb, ':', label = 'Chebyshev order 10')
plt.legend(loc = 'upper right')
plt.xlim(xmin = -1, xmax = 1)
plt.title('Chebyshev polynomial approx. - Runge Function', size=10)
plt.ylabel('f(x)', size = 9)
plt.xlabel('x', size = 9)
plt.show()

plt.plot(roots, y3, label = 'Original function')
plt.plot(roots, a33cheb,'-', label = 'Chebyshev order 3')
plt.plot(roots, a35cheb, '--', label = 'Chebyshev order 5')
plt.plot(roots, a310cheb, ':', label = 'Chebyshev order 10')
plt.legend(loc = 'upper right')
plt.xlim(xmin = -1, xmax = 1)
plt.title('Chebyshev polynomial approx. - Ramp Function', size=10)
plt.ylabel('f(x)', size = 9)
plt.xlabel('x', size = 9)
plt.show()

#We compute the associated errors terms Chebyshev interpolation with Chebyshev polynomial:
e13cheb = abs(y1-a13cheb)
e23cheb = abs(y2-a23cheb)
e33cheb = abs(y3_eval-a33cheb)
e15cheb = abs(y1_eval-a15cheb)
e25cheb = abs(y2_eval-a25cheb)
e35cheb = abs(y3_eval-a35cheb)
e110cheb = abs(y1_eval-a110cheb)
e210cheb = abs(y2_eval-a210cheb)
e310cheb = abs(y3_eval-a310cheb)

#Now, we do the plot for the error terms for the first function:
plt.plot(roots, e13cheb,'-', label = 'Error of order 3')
plt.plot(roots, e15cheb,'--', label = 'Error of order 5')
plt.plot(roots, e110cheb, ':', label = 'Error of order 10')
plt.legend(loc = 'upper right')
plt.xlim(xmin = -1, xmax = 1)
plt.title('Errors Chebyshev approximation with Chebyshev polynomial - Function 1', size=10)
plt.ylabel('f(x)', size=9)
plt.xlabel('x', size=9)
plt.show()

#Now, we do the plot for the error terms for the runge function:
plt.plot(roots, e23cheb,'-', label = 'Error of order 3')
plt.plot(roots, e25cheb,'--', label = 'Error of order 5')
plt.plot(roots, e210cheb, ':', label = 'Error of order 10')
plt.legend(loc = 'upper right')
plt.xlim(xmin = -1, xmax = 1)
plt.title('Errors Chebyshev approximation with Chebyshev polynomial - Runge Function', size=10)
plt.ylabel('f(x)', size=9)
plt.xlabel('x', size=9)
plt.show()

#Now, we do the plot for the error terms for the ramp function:
plt.plot(roots, e33cheb,'-', label = 'Error of order 3')
plt.plot(roots, e35cheb,'--', label = 'Error of order 5')
plt.plot(roots, e310cheb, ':', label = 'Error of order 10')
plt.legend(loc = 'upper right')
plt.xlim(xmin = -1, xmax = 1)
plt.title('Errors Chebyshev approximation with Chebyshev polynomial - Ramp Function', size=10)
plt.ylabel('f(x)', size=9)
plt.xlabel('x', size=9)
plt.show()

