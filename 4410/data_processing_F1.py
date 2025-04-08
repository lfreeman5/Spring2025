import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# --- Previously defined functions ---
def voltage_to_force(lift_voltage, drag_voltage):
    """Convert load cell voltages to forces in Newtons."""
    F_L = 53.411 * lift_voltage - 0.2353
    F_D = 75.424 * drag_voltage + 0.0971
    return F_L, F_D

def compute_coefficients(V, F_L, F_D, S, rho=1.0):
    """Compute lift and drag coefficients."""
    q = 0.5 * rho * V**2  # Dynamic pressure
    C_L = F_L / (q * S)
    C_D = F_D / (q * S)
    return C_L, C_D

def reynolds_number(V, L=0.1, nu=1.49e-5):
    return (V * L) / nu

def plot_coeff_vs_reynolds(data, S, lift=True):
    """Plot Lift Coefficient (C_L) or Drag Coefficient (C_D) as a function of Reynolds Number (Re)."""
    plt.figure()
    for i, (H, dataset) in enumerate(data.items()):
        V, L_volt, D_volt = dataset.T
        F_L, F_D = voltage_to_force(L_volt, D_volt)
        C_L, C_D = compute_coefficients(V, F_L, F_D,S)
        Re = reynolds_number(V)
        color = grayscale_styles[i % len(grayscale_styles)]
        marker = markers[i % len(markers)]
        if lift:
            plt.plot(Re, C_L, color, marker=marker, label=f"H = {H:.1f} mm", markersize=5)
            plt.ylabel("Lift Coefficient, C_L")
            plt.title("Lift Coefficient vs. Reynolds Number")
        else:
            plt.plot(Re, C_D, color, marker=marker, label=f"H = {H:.1f} mm", markersize=5)
            plt.ylabel("Drag Coefficient, C_D")
            plt.title("Drag Coefficient vs. Reynolds Number")
    plt.xlabel("Reynolds Number, Re")
    plt.legend()
    plt.xscale("log")
    plt.tight_layout()
    plt.show()

def calculate_drag_and_velocity(C_D_wing, A_wing):
    """
    Calculate total drag force (F_D) and maximum straight-line velocity (V_s).
    """
    total_C_D = C_D_car * A_car + C_D_wing * A_wing
    # F_D = 0.5 * rho * P / (total_C_D)
    F_D = 0.5 * rho * P / total_C_D  # as given in the original formulation
    V_s = np.cbrt(2 * P / (rho * total_C_D))
    return F_D, V_s

def calculate_uncertainty_drag_velocity(C_D_wing, A_wing, sigma_CD_wing):
    """
    Calculate uncertainty in total drag force (F_D) and maximum straight-line velocity (V_s).
    """
    total_C_D = C_D_car * A_car + C_D_wing * A_wing
    # Using: sigma_{F_D} = (0.5*rho*P*A_wing/(total_C_D^2))*sigma_CD_wing
    sigma_F_D = 0.5 * rho * P * A_wing / (total_C_D**2) * sigma_CD_wing
    # Using: sigma_{V_s} = (A_wing/3)*(2P/rho)**(1/3)*(total_C_D)**(-4/3)*sigma_CD_wing
    sigma_V_s = (A_wing / 3) * (2 * P / rho)**(1/3) * (total_C_D)**(-4/3) * sigma_CD_wing
    return sigma_F_D, sigma_V_s

def calculate_cornering(a, m, g, b, C_L_car, C_L_wing, A_car, A_wing, r, rho):
    """
    Calculate the cornering force and cornering velocity using the updated equation:
    
      (m V_c^2)/r = a*(0.5*rho*V_c^2*(C_L_car*A_car + C_L_wing*A_wing) + m*g) + b
    
    Solving for V_c gives:
    
      V_c = sqrt((a*m*g + b) / ((m/r) - (a*rho*(C_L_car*A_car + C_L_wing*A_wing))/2))
    
    The cornering force is defined as:
    
      F_corner = m*V_c^2 / r.
    """
    C_sum = C_L_car * A_car + C_L_wing * A_wing
    denominator = (m / r) - (a * rho * C_sum) / 2
    numerator = a * m * g + b
    V_c = np.sqrt(numerator / denominator)
    F_corner = m * V_c**2 / r
    return F_corner, V_c

def calculate_uncertainty_cornering(a, m, g, b, C_L_car, C_L_wing, A_car, A_wing, r, rho, sigma_C_L_wing):
    """
    Calculate the uncertainty in cornering force and cornering velocity using the updated equations,
    considering only the uncertainty in C_L_wing.
    
    The cornering velocity is given by:
    
      V_c = sqrt((a*m*g + b) / ((m/r) - (a*rho*(C_L_car*A_car + C_L_wing*A_wing))/2))
    
    and its derivative with respect to C_L_wing is:
    
      dV_c/d(C_L_wing) = (a*rho*A_wing/(4)) * sqrt(a*m*g+b) / ((m/r - (a*rho*(C_L_car*A_car+C_L_wing*A_wing))/2)^(3/2)).
    
    Then, the uncertainty is:
    
      sigma_V_c = |dV_c/d(C_L_wing)| * sigma_C_L_wing,
    
    and for the cornering force:
    
      sigma_F_corner = (2*m*V_c/r)*sigma_V_c.
    """
    C_sum = C_L_car * A_car + C_L_wing * A_wing
    denominator = (m / r) - (a * rho * C_sum) / 2
    numerator = a * m * g + b
    V_c = np.sqrt(numerator / denominator)
    # Compute derivative dV_c/d(C_L_wing)
    dV_c_dCLwing = (a * rho * A_wing / 4) * np.sqrt(numerator) / (denominator**(3/2))
    sigma_V_c = abs(dV_c_dCLwing) * sigma_C_L_wing
    sigma_F_corner = (2 * m * V_c / r) * sigma_V_c
    return sigma_F_corner, sigma_V_c

def compute_uncertainty(V, L_volt, D_volt, S, sigma_Voltage, sigma_V):
    """Calculate the uncertainty in C_L and C_D."""
    F_L = 53.411 * L_volt - 0.2353
    F_D = 75.424 * D_volt + 0.0971
    C_L = (2 * F_L) / (rho * S * V**2)
    C_D = (2 * F_D) / (rho * S * V**2)
    dCL_dV = -2 * F_L / (rho * S * V**3)
    dCL_dLvolt = (2 * 53.411) / (rho * S * V**2)
    dCD_dV = -2 * F_D / (rho * S * V**3)
    dCD_dDvolt = (2 * 75.424) / (rho * S * V**2)
    sigma_CL = np.sqrt((dCL_dV * sigma_V)**2 + (dCL_dLvolt * sigma_Voltage)**2)
    sigma_CD = np.sqrt((dCD_dV * sigma_V)**2 + (dCD_dDvolt * sigma_Voltage)**2)
    return sigma_CL, sigma_CD

def calc_laptime(v_s, v_c, ss, sc):
    sigma = np.sqrt((800/v_s/v_s * ss)**2 + (800/v_c/v_c * sc)**2)
    return 800 * (1 / v_s + 1 / v_c), sigma

# --- Experimental data ---
data = {
    1.2: np.array([
        [4.533531077, 0.002757164, 0.00686729],
        [9.599966542, -0.000409627, 0.021426896],
        [14.51460046, -4.12615E-05, 0.048728087],
        [19.42286432, 0.019331361, 0.096138797],
        [24.34477403, 0.02764709, 0.130923751],
        [29.21348022, 0.099747057, 0.202170073]
    ]),
    19.6: np.array([
        [4.548977, 0.012852085, 0.003833619],
        [9.461271819, 0.046031528, 0.019473319],
        [14.36827689, 0.099366798, 0.042925705],
        [19.30348871, 0.170325322, 0.074424288],
        [24.22476961, 0.267020184, 0.116348812],
        [29.09881407, 0.382952876, 0.16685801]
    ]),
    40.4: np.array([
        [4.474372389, 0.010295755, 0.00452288],
        [9.445136276, 0.038123914, 0.019732452],
        [14.38478578, 0.080196384, 0.043766806],
        [19.27672113, 0.133867338, 0.079635114],
        [24.1733239, 0.211501995, 0.12297035],
        [29.08692583, 0.306117771, 0.172612893]
    ]),
    60.0: np.array([
        [4.547728105, 0.008457657, 0.003493671],
        [9.420518455, 0.028869062, 0.020298246],
        [14.31825203, 0.059706818, 0.044794939],
        [19.19934035, 0.111810639, 0.079287635],
        [24.09324321, 0.172900963, 0.124283959],
        [28.99450315, 0.2722373, 0.178945052]
    ]),
    80.2: np.array([
        [4.587835378, 0.008340249, 0.001488833],
        [9.37046756, 0.024916546, 0.019439821],
        [14.36138077, 0.054993578, 0.044261047],
        [19.18751975, 0.093964314, 0.077856322],
        [24.04702508, 0.157024524, 0.127030333],
        [28.93730143, 0.235267217, 0.186743269]
    ]),
    100.2: np.array([
        [4.572262199, 0.008580978, 0.003093645],
        [9.331366237, 0.023463414, 0.018720293],
        [14.28246922, 0.04794024, 0.046377122],
        [19.12485481, 0.0846846, 0.080294357],
        [24.04511754, 0.145879654, 0.120838853],
        [28.96612215, 0.218736898, 0.175103716]
    ])
}

# --- Matplotlib settings ---
mpl.rcParams.update({
    'font.family': 'Times New Roman',
    'font.size': 10,
    'axes.labelsize': 10,
    'axes.titlesize': 12,
    'legend.fontsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'lines.linewidth': 1,
    'lines.markersize': 4,
    'grid.linestyle': '--',
    'grid.linewidth': 0.5,
    'axes.grid': True,
    'axes.grid.which': 'both',
    'axes.linewidth': 0.75,
    'xtick.major.size': 5,
    'ytick.major.size': 5,
    'xtick.minor.size': 2.5,
    'ytick.minor.size': 2.5,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'legend.loc': 'best',
    'legend.frameon': False,
    'savefig.dpi': 300,
    'savefig.format': 'pdf',
    'figure.figsize': (6.5, 4.5),
    'figure.dpi': 100
})
grayscale_styles = ['black', 'dimgray', 'gray', 'darkgray']
markers = ['o', 's', '^', 'd', 'v', '*']

if __name__ == "__main__":
    # Constants
    rho = 1.0  # Air density (kg/m³)
    A_wing_tunnel = 0.06693  # Wing area from tunnel test (m²)
    A_wing = 9 * A_wing_tunnel  # Scaled wing area (m²)
    A_car = 3.0  # Car frontal area (m²)
    C_D_car = 1.0  # Drag coefficient of car
    C_L_car = 0.5  # Lift coefficient of car
    P = 74569.99  # Power in watts
    m = 200  # Mass of vehicle (kg)
    g = 9.80  # Gravity (m/s²)
    r = 127.0  # Turn radius (m)
    a = 0.198  # Cornering force constant
    b = 1      # Cornering force offset
    sigma_V = 0.1  # Velocity uncertainty (m/s)
    sigma_Voltage = 0.02  # Voltage uncertainty (V)

    for H, dataset in data.items():
        V, L_volt, D_volt = dataset.T
        F_L, F_D = voltage_to_force(L_volt, D_volt)
        C_L, C_D = compute_coefficients(V, F_L, F_D,A_wing_tunnel)

        # Use highest velocity at each height
        idx_max_V = np.argmax(V)
        V_max = V[idx_max_V]
        C_L_max = C_L[idx_max_V]
        C_D_max = C_D[idx_max_V]
        sigma_C_L_wing, sigma_C_D_wing = compute_uncertainty(V_max, L_volt[idx_max_V], D_volt[idx_max_V], A_wing_tunnel, sigma_Voltage, sigma_V)
        print(f'Uncertainty in C_L, C_D: {sigma_C_L_wing} {sigma_C_D_wing}')

        # Compute drag force and straight-line velocity
        F_D_total, V_s = calculate_drag_and_velocity(C_D_max, A_wing)
        
        # Compute cornering force and cornering velocity
        F_corner, V_c = calculate_cornering(a, m, g, b, C_L_car, C_L_max, A_car, A_wing, r, rho)

        # Compute uncertainties
        sigma_F_D_total, sigma_V_s = calculate_uncertainty_drag_velocity(C_D_max, A_wing, sigma_C_D_wing)
        sigma_F_corner, sigma_V_c = calculate_uncertainty_cornering(a, m, g, b, C_L_car, C_L_max, A_car, A_wing, r, rho, sigma_C_L_wing)

        # Print results
        print(f"Height: {H:.1f} mm")
        print(f"  Max Velocity: {V_max:.2f} m/s")
        print(f"  C_L: {C_L_max:.4f}, C_D: {C_D_max:.4f}")
        print(f"  Total Drag Force: {F_D_total:.4f} ± {sigma_F_D_total:.4f}")
        print(f"  Straight-Line Velocity: {V_s:.4f} ± {sigma_V_s:.4f}")
        print(f"  Cornering Force: {F_corner:.4f} ± {sigma_F_corner:.4f}")
        print(f"  Cornering Velocity: {V_c:.4f} ± {sigma_V_c:.4f}")
        l, sl = calc_laptime(V_s, V_c, sigma_V_s, sigma_V_c)
        print(f"  Lap Time: {l:.4f} ± {sl:.4f}")
        print()

    # Lap Time without wing:
    _, V_s_no_wing = calculate_drag_and_velocity(0,0)
    _, V_c_no_wing = calculate_cornering(a, m, g, b, 0, C_L_max, A_car, 0, r, rho)
    t_no_wing = calc_laptime(V_s_no_wing,V_c_no_wing,0,0)
    print(f'Straight line & cornering vel w/o wing: {V_s_no_wing}     {V_c_no_wing}')
    print(f'Lap time with no wing: {t_no_wing}')

    # # Generate plots
    plot_coeff_vs_reynolds(data, A_wing_tunnel)
    plot_coeff_vs_reynolds(data, A_wing_tunnel, lift=False)

    heights = [1.2, 19.6, 40.4, 60.0, 80.2, 100.2]  # Heights in mm
    speeds = [5, 10, 15, 20, 25, 30]  # Set speeds in m/s

    def print_table(force_matrix, force_name):
        print(f"{force_name} (N)")
        print("Speed (m/s)\t" + "\t".join([f"H={h} mm" for h in heights]))
        for i, V in enumerate(speeds):
            row = [f"{V}"] + [f"{force_matrix[i, j]:.2f}" for j in range(len(heights))]
            print("\t".join(row))
        print("\n")

    # Compute forces
    L_matrix = np.zeros((len(speeds), len(heights)))
    D_matrix = np.zeros((len(speeds), len(heights)))

    for j, H in enumerate(heights):
        L, D = voltage_to_force(*data[H][:, 1:].T)  # Convert voltages to forces
        for i, V in enumerate(speeds):
            idx = np.argmin(np.abs(data[H][:, 0] - V))  # Find closest velocity
            L_matrix[i, j] = L[idx]
            D_matrix[i, j] = D[idx]

    # Print tables
    print_table(L_matrix, "Lift Force")
    print_table(D_matrix, "Drag Force")

    # Given parameters
    heights = np.array([1.2, 19.6, 40.4, 60.0, 80.2, 100.2])  # Heights in mm
    c = 10  # Reference chord length in mm
    H_non_dim = heights / c  # Non-dimensionalized height
    speeds = np.array([5, 10, 15, 20, 25, 30])  # Set speeds in m/s

    # Example force data (replace with actual data)
    S_ref = 0.067  # Reference area (m^2)

    C_L_matrix = np.zeros((6, 6))  # Placeholder for computed lift coefficients
    for j, height in enumerate(heights):  # Loop over heights
        for i in range(6):  # Loop over velocities
            V, F_L, F_D = data[height][i]
            _, C_L_matrix[i, j] = compute_coefficients(V, F_L, F_D, S_ref)

    # Grayscale styles and markers
    grayscale_styles = ['black', 'dimgray', 'gray', 'darkgray', 'silver', 'lightgray']
    markers = ['o', 's', '^', 'd', 'v', '*']

    # Plot
    plt.figure(figsize=(8, 6))
    for i, V in enumerate(speeds):
        plt.plot(H_non_dim, C_L_matrix[i, :], 
                marker=markers[i], color=grayscale_styles[i], 
                label=f"V = {V} m/s", linestyle='-')

    plt.xlabel("Non-dimensional Height (H/c)")
    plt.ylabel("Drag Coefficient $C_D$")
    plt.title("Drag Coefficient vs. Non-dimensional Height")
    plt.legend()
    plt.grid(True)
    plt.show()

    wing_heights = np.array([1.2, 19.6, 40.4, 60.0, 80.2, 100.2])
    lap_times = np.array([70.7689, 70.0991, 70.2688, 70.3583, 70.4596, 70.4485])
    errors = np.array([0.0813, 0.0829, 0.0828, 0.0831, 0.0832, 0.0833])

    plt.figure(figsize=(8,6))
    plt.errorbar(wing_heights, lap_times, yerr=errors, fmt='-o', color='black', ecolor='black', capsize=5)
    plt.xlabel("Wing height (mm)")
    plt.ylabel("Lap Time (s)")
    plt.show()