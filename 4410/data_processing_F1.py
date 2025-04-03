import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def voltage_to_force(lift_voltage, drag_voltage):
    """Convert load cell voltages to forces in Newtons."""
    F_L = 53.411 * lift_voltage - 0.2353
    F_D = 75.424 * drag_voltage + 0.0971
    return F_L, F_D

def compute_coefficients(V, F_L, F_D, rho=1.0, S=0.1):
    """Compute lift and drag coefficients."""
    q = 0.5 * rho * V**2  # Dynamic pressure
    C_L = F_L / (q * S)
    C_D = F_D / (q * S)
    return C_L, C_D

def reynolds_number(V, L=0.1, nu=1.49e-5):
    return (V * L) / nu

def plot_coeff_vs_reynolds(data, lift=True):
    """Plot Lift Coefficient (C_L) or Drag Coefficient (C_D) as a function of Reynolds Number (Re)."""
    plt.figure()

    for i, (H, dataset) in enumerate(data.items()):
        V, L_volt, D_volt = dataset.T
        F_L, F_D = voltage_to_force(L_volt, D_volt)
        C_L, C_D = compute_coefficients(V, F_L, F_D)
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
    plt.xscale("log")  # Reynolds number is often plotted on a log scale
    plt.tight_layout()
    plt.show()

sigma_V = 0.1  # Uncertainty in velocity (m/s)
sigma_Voltage = 0.02  # Uncertainty in load cell voltage (V)
rho = 1.0  # Air density (kg/m^3)
S = 0.06693  # Wing area (m^2)

P = 74569.99  # Power in watts (example value, adjust as needed)

# Example constants for the car (should be provided)
C_D_car = 1  # Drag coefficient of car (example value)
A_car = 3.0  # Cross-sectional area of car (m^2)

def calculate_drag_and_velocity(C_D_wing, A_wing):
    """
    Calculate total drag force (F_D) and maximum straight-line velocity (V_s).
    
    Parameters:
    C_D_wing (float): Drag coefficient of the wing
    A_wing (float): Area of the wing (m^2)
    
    Returns:
    F_D (float): Total drag force in newtons (N)
    V_s (float): Maximum straight-line velocity in m/s
    """
    # Calculate total drag coefficient
    total_C_D = C_D_car * A_car + C_D_wing * A_wing
    
    # Calculate drag force
    # Assume the vehicle is moving at V_s, so calculate drag force at V_s
    F_D = 0.5 * rho * P / (total_C_D)  # Using the equation for total drag force
    
    # Calculate maximum straight-line velocity
    V_s = np.cbrt(2 * P / (rho * total_C_D))  # Using the equation for V_s
    
    return F_D, V_s

def calculate_uncertainty_drag_velocity(C_D_wing, A_wing, sigma_CD_wing):
    """
    Calculate uncertainty in total drag force (F_D) and maximum straight-line velocity (V_s).
    
    Parameters:
    C_D_wing (float): Drag coefficient of the wing
    A_wing (float): Area of the wing (m^2)
    sigma_CD_wing (float): Uncertainty in the drag coefficient of the wing
    
    Returns:
    sigma_F_D (float): Uncertainty in total drag force (N)
    sigma_V_s (float): Uncertainty in maximum velocity (m/s)
    """
    # Calculate total drag coefficient (drag coefficients of car + wing)
    total_C_D = C_D_car * A_car + C_D_wing * A_wing
    
    # Calculate maximum straight-line velocity
    V_s = np.cbrt(2 * P / (rho * total_C_D))
    
    # Uncertainty in drag force (F_D)
    sigma_F_D = 0.5 * rho * V_s**2 * A_wing * sigma_CD_wing
    
    # Uncertainty in velocity (V_s)
    sigma_V_s = (2 * P / (3 * rho * total_C_D**(4/3))) * sigma_CD_wing
    
    return sigma_F_D, sigma_V_s

def calculate_cornering(a, m, g, b, C_L_car, C_L_wing, A_car, A_wing, r, rho):
    """
    Calculate the cornering force and cornering velocity.
    
    Parameters:
    a (float): Constant 'a'
    m (float): Mass of the car (kg)
    g (float): Gravitational acceleration (m/s^2)
    b (float): Constant 'b'
    C_L_car (float): Lift coefficient of the car
    C_L_wing (float): Lift coefficient of the wing
    A_car (float): Area of the car (m^2)
    A_wing (float): Area of the wing (m^2)
    r (float): Radius of the curve (m)
    rho (float): Air density (kg/m^3)
    
    Returns:
    F_corner (float): Cornering force (N)
    V_c (float): Cornering velocity (m/s)
    """
    
    # Cornering velocity calculation
    denominator = (m / r) + (a / 2) * rho * (C_L_car * A_car + C_L_wing * A_wing)
    numerator = a * m * g + b
    V_c = np.sqrt(numerator / denominator)
    
    # Cornering force calculation
    F_corner = a * (-0.5 * rho * V_c**2 * (C_L_car * A_car + C_L_wing * A_wing) + m * g) + b
    
    return F_corner, V_c

def calculate_uncertainty_cornering(a, m, g, b, C_L_car, C_L_wing, A_car, A_wing, r, rho, sigma_C_L_wing):
    """
    Calculate the uncertainty in cornering force and cornering velocity, considering only the uncertainty in C_L_wing.
    
    Parameters:
    sigma_C_L_wing (float): Uncertainty in C_L_wing.
    (other parameters as before)

    Returns:
    sigma_F_corner (float): Uncertainty in cornering force (N)
    sigma_V_c (float): Uncertainty in cornering velocity (m/s)
    """
    
    # Calculate cornering velocity V_c (as before)
    denominator = (m / r) + (a / 2) * rho * (C_L_car * A_car + C_L_wing * A_wing)
    numerator = a * m * g + b
    V_c = np.sqrt(numerator / denominator)
    
    # Calculate the partial derivative of cornering velocity with respect to C_L_wing
    partial_V_c_C_L_wing = - (a * (m * g + b)) / (2 * r * (denominator**2)) * (A_wing * rho)  # derivative with respect to C_L_wing
    
    # Uncertainty in cornering velocity
    sigma_V_c = abs(partial_V_c_C_L_wing) * sigma_C_L_wing
    
    # Now calculate the partial derivative of cornering force with respect to C_L_wing
    partial_F_corner_C_L_wing = -0.5 * rho * V_c**2 * A_wing * a  # derivative with respect to C_L_wing
    
    # Uncertainty in cornering force
    sigma_F_corner = abs(partial_F_corner_C_L_wing) * sigma_C_L_wing
    
    return sigma_F_corner, sigma_V_c


def compute_uncertainty(V, L_volt, D_volt):
    """Calculate the uncertainty in C_L and C_D."""
    
    # Compute forces from voltages
    F_L = 53.411 * L_volt - 0.2353
    F_D = 75.424 * D_volt + 0.0971
    
    # Compute C_L and C_D
    C_L = (2 * F_L) / (rho * S * V**2)
    C_D = (2 * F_D) / (rho * S * V**2)
    
    # Compute partial derivatives for uncertainty propagation
    dCL_dV = -2 * F_L / (rho * S * V**3)
    dCL_dLvolt = (2 * 53.411) / (rho * S * V**2)

    dCD_dV = -2 * F_D / (rho * S * V**3)
    dCD_dDvolt = (2 * 75.424) / (rho * S * V**2)

    # Apply uncertainty propagation formula
    sigma_CL = np.sqrt((dCL_dV * sigma_V)**2 + (dCL_dLvolt * sigma_Voltage)**2)
    sigma_CD = np.sqrt((dCD_dV * sigma_V)**2 + (dCD_dDvolt * sigma_Voltage)**2)

    return sigma_CL, sigma_CD

# Experimental data grouped by height (H, V, L, D)
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

mpl.rcParams.update({
    'font.family': 'Times New Roman',  # AIAA recommends Times New Roman font
    'font.size': 10,                   # Standard font size for AIAA documents
    'axes.labelsize': 10,              # Font size for axis labels
    'axes.titlesize': 12,              # Font size for plot titles
    'legend.fontsize': 8,              # Font size for legends
    'xtick.labelsize': 8,              # Font size for x-axis tick labels
    'ytick.labelsize': 8,              # Font size for y-axis tick labels
    'lines.linewidth': 1,              # Line width for plot lines
    'lines.markersize': 4,             # Marker size for plot markers
    'grid.linestyle': '--',            # Dashed grid lines
    'grid.linewidth': 0.5,             # Grid line width
    'axes.grid': True,                 # Enable grid by default
    'axes.grid.which': 'both',         # Apply grid to both major and minor ticks
    'axes.linewidth': 0.75,            # Axis line width
    'xtick.major.size': 5,             # Major tick size for x-axis
    'ytick.major.size': 5,             # Major tick size for y-axis
    'xtick.minor.size': 2.5,           # Minor tick size for x-axis
    'ytick.minor.size': 2.5,           # Minor tick size for y-axis
    'xtick.direction': 'in',           # Inward ticks for x-axis
    'ytick.direction': 'in',           # Inward ticks for y-axis
    'legend.loc': 'best',              # Optimal legend placement
    'legend.frameon': False,           # No frame around legend
    'savefig.dpi': 300,                # High resolution for saved figures
    'savefig.format': 'pdf',           # Save figures in PDF format
    'figure.figsize': (6.5, 4.5),      # Figure size in inches (width, height)
    'figure.dpi': 100                  # Display resolution
})

grayscale_styles = ['black', 'dimgray', 'gray', 'darkgray']
markers = ['o', 's', '^', 'd', 'v', '*']

if __name__ == "__main__":
    # Constants
    rho = 1.0  # Air density (kg/m³)
    A_wing = 0.06693  # Wing area (m²)
    A_car = 2.0  # Car frontal area (m²)
    C_D_car = 0.3  # Drag coefficient of car
    C_L_car = 0.1  # Lift coefficient of car
    P = 74569.99  # Power in watts
    m = 200  # Mass of vehicle (kg)
    g = 9.80  # Gravity (m/s²)
    r = 127.0  # Turn radius (m)
    a = 0.198  # Cornering force constant
    b = 1  # Cornering force offset


    for H, dataset in data.items():
        V, L_volt, D_volt = dataset.T
        F_L, F_D = voltage_to_force(L_volt, D_volt)
        C_L, C_D = compute_coefficients(V, F_L, F_D)

        # Use highest velocity at each height
        idx_max_V = np.argmax(V)
        V_max = V[idx_max_V]
        C_L_max = C_L[idx_max_V]
        C_D_max = C_D[idx_max_V]
        sigma_C_L_wing, sigma_C_D_wing = compute_uncertainty(V_max,L_volt[idx_max_V],D_volt[idx_max_V])
        print(f'Uncertainty in CL, CD: {sigma_C_L_wing} {sigma_C_D_wing}')



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
        print(f"  Total Drag Force: {F_D_total:.4f} N ± {sigma_F_D_total:.4f} N")
        print(f"  Straight-Line Velocity: {V_s:.4f} m/s ± {sigma_V_s:.4f} m/s")
        print(f"  Cornering Force: {F_corner:.4f} N ± {sigma_F_corner:.4f} N")
        print(f"  Cornering Velocity: {V_c:.4f} m/s ± {sigma_V_c:.4f} m/s")
        print()

    # Generate plots
    plot_coeff_vs_reynolds(data)
    plot_coeff_vs_reynolds(data, lift=False)
