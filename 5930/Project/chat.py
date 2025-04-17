import cantera as ct

def compute_laminar_flame_speed(phi, p_atm, T0,
                                mech='./chem_files/nh3_san_diego.yaml'):
    """
    Compute the laminar flame speed for NH3/air using Cantera.

    Parameters
    ----------
    phi : float
        Equivalence ratio (fuel-to-oxidizer ratio relative to stoichiometric).
    p_atm : float
        Pressure in atmospheres.
    T0 : float
        Unburned gas temperature in kelvin.
    mech : str, optional
        Path to the Cantera mechanism YAML file (default './chem_files/nh3_san_diego.yaml').

    Returns
    -------
    float
        Laminar flame speed in m/s.
    """
    # 1. Create gas mixture
    gas = ct.Solution(mech)  # Load mechanism :contentReference[oaicite:10]{index=10}
    gas.set_equivalence_ratio(phi, 'NH3', {'O2':1, 'N2':3.76})  # Set φ :contentReference[oaicite:11]{index=11}
    gas.TP = T0, p_atm * ct.one_atm  # Set initial state :contentReference[oaicite:12]{index=12}

    # 2. Set up 1-D freely-propagating flame
    flame = ct.FreeFlame(gas, width=0.03)  # Flat flame domain :contentReference[oaicite:13]{index=13}
    flame.set_refine_criteria(ratio=3.0, slope=0.06, curve=0.12)  # Grid refinement :contentReference[oaicite:14]{index=14}

    # 3. Solve for steady flame structure
    flame.solve(loglevel=0, auto=True)  # Solve flame :contentReference[oaicite:15]{index=15}

    # 4. Extract laminar flame speed (inlet velocity)
    return flame.velocity[0]  # Laminar flame speed [m/s] :contentReference[oaicite:16]{index=16}

if __name__ == "__main__":
    Su = compute_laminar_flame_speed(phi=1.0, p_atm=1.0, T0=300.0)
    print(f"Laminar flame speed at φ=1.0, p=1 atm, T₀=300 K: {Su:.3f} m/s")
