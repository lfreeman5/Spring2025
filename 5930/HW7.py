import cantera as ct

# Standard state conditions
T = 298.15  # K
P = ct.one_atm  # Pa

# Load n-hexane from 'chemicals.yaml'
n_hexane = ct.Solution('chemicals.yaml', 'n-hexane')
n_hexane.TP = T, P

# Load air from 'air.yaml'
air = ct.Solution('air.yaml')
air.TP = T, P

# Get properties
hf_nhexane = n_hexane.enthalpy_mole          # J/mol
cp_nhexane = n_hexane.cp_mass                # J/kg-K
k_nhexane = n_hexane.thermal_conductivity    # W/m-K

k_air = air.thermal_conductivity             # W/m-K

# Output
print(f"n-Hexane at {T} K and 1 atm:")
print(f"  Enthalpy of formation (h_f): {hf_nhexane:.2f} J/mol")
print(f"  Specific heat (cp): {cp_nhexane:.2f} J/kg-K")
print(f"  Thermal conductivity (k): {k_nhexane:.5f} W/m-K\n")

print(f"Air at {T} K and 1 atm:")
print(f"  Thermal conductivity (k): {k_air:.5f} W/m-K")
