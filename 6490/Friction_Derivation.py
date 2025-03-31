import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib.pyplot as plt

nu = 1e-4
kappa = 0.41
h = 1
u_cl = 10
rho=1
mu = rho*nu

def ode_system(x, u, tau_w):
    du_dx = (-1 + np.sqrt(1 + 4 * (tau_w / (nu*mu)) * kappa**2 * x**2)) / (2 * kappa**2 * x**2 / nu)
    return du_dx

def objective(tau_w):
    def ode(t, u):
        return ode_system(t, u, tau_w)
    
    initial_condition = [0.0]
    solution = solve_ivp(ode, (1e-7, h), initial_condition, t_eval=np.linspace(1e-7, h, 1000))
    
    u_at_x_target = solution.y[0][-1]
    return u_at_x_target - u_cl

result = brentq(objective, 1e-5, 1e+6)

# After finding tau_w, compute the solution over the entire interval and calculate the average u
def compute_average_u(tau_w):
    def ode(t, u):
        return ode_system(t, u, tau_w)
    
    initial_condition = [0.0]
    solution = solve_ivp(ode, (1e-7, h), initial_condition, t_eval=np.linspace(1e-7, h, 1000))
    
    u_values = solution.y[0]
    x_values = solution.t
    average_u = np.trapz(u_values, x_values) / (x_values[-1] - x_values[0])  # Trapezoidal integration
    return average_u

average_u = compute_average_u(result)

print(f"Found tau_w: {result}")
print(f"Average value of u(x) over the interval: {average_u}")

cf_theory = 2*result / (rho*u_cl**2)
Re_b = 2*h*average_u/nu
cf_exp = 0.073*(Re_b)**-0.25

print(f'Theory value for cf: {cf_theory}')
print(f'Experimental value for cf: {cf_exp}')