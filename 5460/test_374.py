import numpy as np
from bifurcation import *

def xdot(x,a,h):
    return x*(1-x)-h*x/(a+x)

A=2
def f(x,h):
    return xdot(x,A,h)

d = bifurcation_search((-0.25,.25),(1.5,2.5),f)
s = characterize_stability(f,d)
plot_bifurcation_data(d,s)    

A=0.5
def f(x,h):
    return xdot(x,A,h)

d = bifurcation_search((-0.1,0.1),(0.4,0.8),f)
s = characterize_stability(f,d)
plot_bifurcation_data(d,s)    