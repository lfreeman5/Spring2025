import cantera as ct
import numpy as np
import copy
import matplotlib.pyplot as plt
import pandas as pd

def combustion_mixture(phi):
    """computes reaction coefficients for nh3, dme, air"""

    nh3_frac = 1
    o2_stoich_per_mol = nh3_frac * 0.75

    o2_actual = o2_stoich_per_mol / phi
    n2_actual = o2_actual*3.76

    return nh3_frac, o2_actual, n2_actual

def plot_flame_structure(flame, gas):
    z = flame.grid  # flame coordinate
    T = flame.T  # temperature profile
    HRR = flame.heat_release_rate  # W/m^3

    species = ['NH3', 'O2', 'N2', 'NO', 'H2O']
    Y = {sp: flame.X[gas.species_index(sp), :] for sp in species}

    plt.figure(figsize=(10, 6))
    plt.subplot(3, 1, 1)
    plt.plot(z, T, label='Temperature (K)', color='red')
    plt.ylabel('Temperature [K]')
    plt.grid()
    plt.title('Flame Structure')
    
    plt.subplot(3, 1, 2)
    plt.plot(z, HRR, label='Heat Release Rate', color='orange')
    plt.ylabel('HRR [W/m³]')
    plt.grid()

    plt.subplot(3, 1, 3)
    for sp in species:
        plt.plot(z, Y[sp], label=sp)
    plt.ylabel('Mole Fraction')
    plt.xlabel('Flame Coordinate [m]')
    plt.legend()
    plt.grid()

    plt.tight_layout()
    plt.show()

def calc_flame_thickness(gas, flame, x_target):
    idx = np.argmin(np.abs(flame.grid - x_target))

    T = flame.T[idx]
    Y = flame.Y[:, idx]
    gas.TPY = T, flame.P, Y
    alpha = gas.thermal_conductivity / (gas.density * gas.cp_mass)
    delta = 2 * alpha / flame.velocity[0]
    return delta


def compute_flame_speed(phi, T0, p_atm,
                                mech='./chem_files/nh3_san_diego.yaml'):
    gas = ct.Solution(mech)  # Load mechanism :contentReference[oaicite:10]{index=10}
    gas.set_equivalence_ratio(phi, 'NH3', {'O2':1, 'N2':3.76})  # Set φ :contentReference[oaicite:11]{index=11}
    gas.TP = T0, p_atm * ct.one_atm  # Set initial state :contentReference[oaicite:12]{index=12}

    flame = ct.FreeFlame(gas, width=0.03)  # Flat flame domain :contentReference[oaicite:13]{index=13}
    flame.set_refine_criteria(ratio=3.0, slope=0.06, curve=0.12)  # Grid refinement :contentReference[oaicite:14]{index=14}
    flame.solve(loglevel=0, auto=True)  # Solve flame :contentReference[oaicite:15]{index=15}

    return flame


if __name__ == "__main__":
    gas=ct.Solution('./chem_files/nh3_san_diego.yaml')

    phi_range = np.linspace(0.8,1.2,10)
    T_range = np.linspace(300,500,10)
    p_range = np.linspace(1,10,10)

    phi_results = np.zeros_like(phi_range)
    T_results = np.zeros_like(T_range)
    p_results = np.zeros_like(p_range)

    phi_default=1
    T_default=300
    p_default=1


    flame_default=compute_flame_speed(phi_default,T_default,p_default)

    # print(f'Laminar flame speed default: {flame_default.velocity[0]}')

    # last_flame = flame_default

    # for i,phi in enumerate(phi_range):
    #     print(i)
    #     f = compute_flame_speed(phi,T_default,p_default)
    #     phi_results[i] = f.velocity[0]
    #     print(phi_results[i])


    # for i,p in enumerate(p_range):
    #     print(i*10)
    #     f = compute_flame_speed(phi_default,T_default,p)
    #     p_results[i] = f.velocity[0]

    # for i,T in enumerate(T_range):
    #     print(i*100)
    #     f = compute_flame_speed(phi_default,T,p_default)
    #     T_results[i] = f.velocity[0]

    # plt.figure(figsize=(6,4))
    # plt.plot(phi_range, phi_results, marker='o', color='crimson')
    # plt.xlabel('Equivalence Ratio (φ)')
    # plt.ylabel('Flame Speed (m/s)')
    # plt.title('Laminar Flame Speed vs Equivalence Ratio')
    # plt.grid(True)
    # plt.tight_layout()
    # plt.savefig('flame_speed_vs_phi.png')

    # # Flame speed vs pressure
    # plt.figure(figsize=(6,4))
    # plt.plot(p_range, p_results, marker='s', color='teal')
    # plt.xlabel('Pressure (atm)')
    # plt.ylabel('Flame Speed (m/s)')
    # plt.title('Laminar Flame Speed vs Pressure')
    # plt.grid(True)
    # plt.tight_layout()
    # plt.savefig('flame_speed_vs_pressure.png')

    # # Flame speed vs temperature
    # plt.figure(figsize=(6,4))
    # plt.plot(T_range, T_results, marker='^', color='darkorange')
    # plt.xlabel('Unburned Temperature (K)')
    # plt.ylabel('Flame Speed (m/s)')
    # plt.title('Laminar Flame Speed vs Temperature')
    # plt.grid(True)
    # plt.tight_layout()
    # plt.savefig('flame_speed_vs_temperature.png')

    # plt.show()

    thickness = calc_flame_thickness(gas, flame_default, 0.0215)
    print(f'Flame thickness: {thickness}')
    plot_flame_structure(flame_default,gas)

    sens = pd.DataFrame(index=gas.reaction_equations(), columns=["sensitivity"])

    sens.sensitivity = flame_default.get_flame_speed_reaction_sensitivities()

    sens.head(10)

    sens = sens.iloc[(-sens['sensitivity'].abs()).argsort()]

    fig, ax = plt.subplots()

    sens.head(15).plot.barh(ax=ax, legend=None)

    ax.invert_yaxis()  # put the largest sensitivity on top
    ax.set_title("Sensitivities for NH3-Air combustion")
    ax.set_xlabel(r"Sensitivity: $\frac{\partial\:\ln S_u}{\partial\:\ln k}$")
    ax.grid(axis='x')
    plt.tight_layout()
    plt.show()