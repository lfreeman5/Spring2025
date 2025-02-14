import numpy as np
from bifurcation import *
from scipy.integrate import odeint

if __name__ == "__main__":
    def f(x,r):
        return x+np.tanh(x*r)
    d = bifurcation_search((-10,10),(-3,1),f)
    s = characterize_stability(f,d)
    plot_bifurcation_data(d,s)    

    def f2(x,t): # Dummy variable t used for scipy ode
        return f(x,-2)

    t=np.linspace(0,10,100)
    x01 = -0.5
    x02 = 0.5
    x1 = odeint(f2,x01,t)
    x2 = odeint(f2,x02,t)


    plt.plot(t,x1,label='X0 = -0.5')
    plt.plot(t,x2,label='X0 = 0.5')
    plt.xlabel(f'Time')
    plt.ylabel(f'X(t)')
    plt.title(f'Solutions to dynamical system')
    plt.legend()
    plt.grid(True)
    plt.show()