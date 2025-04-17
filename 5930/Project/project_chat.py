import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

DME_FRACTION = 0.6

def combustion_mixture(phi, dme_frac):
    nh3 = 1 - dme_frac
    o2_stoich = nh3*0.75 + dme_frac*3.0
    air_stoich = o2_stoich*(1+3.76)
    return nh3, dme_frac, air_stoich/phi

# ----------------------------------------------------
# 1) INITIAL SETUP & FIRST SOLVE (with auto=True)
# ----------------------------------------------------
gas = ct.Solution('./chem_files/nh3_dme_mechanism.yaml')
phi0, T0, p0 = 1.0, 300, 1.0
nh3, dme, air = combustion_mixture(phi0, DME_FRACTION)
react0 = f'NH3:{nh3}, CH3OCH3:{dme}, N2:{air*3.76/4.76}, O2:{air/4.76}'
gas.TPX = T0, p0*ct.one_atm, react0

flame = ct.FreeFlame(gas, width=0.05)             # smaller domain
flame.set_refine_criteria(ratio=5, slope=0.1, curve=0.2)
flame.max_grid_points = 200                    # cap grid
flame.transport_model = 'Mix'

# do the only mesh‐refine solve
flame.solve(loglevel=1, auto=True)
print(f"Base speed: {flame.velocity[0]:.4f} m/s")

# ----------------------------------------------------
# 2) PARAMETER SWEEPS (reuse the mesh, auto=False)
# ----------------------------------------------------
phi_range = np.linspace(0.8, 1.2, 10)
p_range   = np.linspace(1,   10, 10)
T_range   = np.linspace(300, 500, 10)

phi_res = np.zeros_like(phi_range)
p_res   = np.zeros_like(p_range)
T_res   = np.zeros_like(T_range)

# Sweep φ
for i, phi in enumerate(phi_range):
    nh3, dme, air = combustion_mixture(phi, DME_FRACTION)
    gas.TPX = T0, p0*ct.one_atm, f'NH3:{nh3}, CH3OCH3:{dme}, N2:{air*3.76/4.76}, O2:{air/4.76}'
    flame.set_initial_guess()
    flame.solve(loglevel=0, auto=False)
    phi_res[i] = flame.velocity[0]

# Sweep pressure
for i, p in enumerate(p_range):
    gas.TPX = T0, p*ct.one_atm, react0
    flame.set_initial_guess()
    flame.solve(loglevel=0, auto=False)
    p_res[i] = flame.velocity[0]

# Sweep temperature
for i, T in enumerate(T_range):
    gas.TPX = T, p0*ct.one_atm, react0
    flame.set_initial_guess()
    flame.solve(loglevel=0, auto=False)
    T_res[i] = flame.velocity[0]

# ----------------------------------------------------
# 3) PLOTTING
# ----------------------------------------------------
plt.figure()
plt.plot(phi_range, phi_res, 'o-')
plt.xlabel('φ'); plt.ylabel('S_L (m/s)')
plt.title('Flame Speed vs Equivalence Ratio')
plt.grid(True)

plt.figure()
plt.plot(p_range, p_res, 's-')
plt.xlabel('Pressure (atm)'); plt.ylabel('S_L (m/s)')
plt.title('Flame Speed vs Pressure')
plt.grid(True)

plt.figure()
plt.plot(T_range, T_res, '^-')
plt.xlabel('Unburned T (K)'); plt.ylabel('S_L (m/s)')
plt.title('Flame Speed vs Temperature')
plt.grid(True)

plt.show()
