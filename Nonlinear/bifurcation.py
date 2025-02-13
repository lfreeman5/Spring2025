import numpy as np
from scipy.optimize import fsolve

def f(x,r):
    return 1-x-np.exp(-r*x)

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def bifurcation_search(x_range, r_range, func):
    (x_min, x_max) = x_range
    (r_min, r_max) = r_range
    
    r_values = np.linspace(r_min, r_max, 500)
    
    bifurcation_data = []
    
    for r in r_values:
        x_stars = set()  # Use a set to store unique x_star values
        def f_r(x):
            return func(x, r)
        
        x_initial = np.linspace(x_min, x_max, 7)
        
        for x_0 in x_initial:
            try:
                x_star = fsolve(f_r, x_0).item()
                # if(r<1):
                #     print(f'x star value: {x_star}')
                if(abs(f_r(x_star))<1e-5):
                    x_stars.add(round(x_star, 6))  # Round to avoid floating-point duplicates
            except:
                pass
        
        bifurcation_data.extend((r, x) for x in x_stars)
    
    bifurcation_data = np.array(bifurcation_data)

    return bifurcation_data
    
def characterize_stability(func, bifurcation_points):
    r_vals = bifurcation_points[:,0]
    x_vals = bifurcation_points[:,1]
    stability_vals = np.zeros_like(x_vals)
    for i,x in enumerate(x_vals):
        stability_vals[i] = cd_estimate(func,x,r_vals[i])
    return stability_vals

def cd_estimate(f,x,r):
    cd_step = 1e-6
    return (f(x+cd_step,r)-f(x-cd_step,r))/(2*cd_step)

def plot_bifurcation_data(data, stability_vals):
    print(data.size)
    if data.size > 0:
        plt.figure(figsize=(8, 6))
        for i in range(len(data[:,0])):
            marker = '.' if stability_vals[i]<0 else 'x'
            color = 'black' if stability_vals[i]<0 else 'red'
            plt.scatter(data[i, 0], data[i, 1], s=4, color = color, marker=marker)
        plt.xlabel("r")
        plt.ylabel("x*")
        plt.title("Bifurcation Diagram")
        plt.show()
    else:
        print("No bifurcation points found.")


if __name__ == "__main__":
    d = bifurcation_search((-10,10),(0,4),f)
    s = characterize_stability(f,d)
    plot_bifurcation_data(d,s)