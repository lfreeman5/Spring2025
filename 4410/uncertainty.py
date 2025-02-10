import numpy as np

def uncertainty_eta(R, cp, Th, Tc, V1, V2, sigma_Th, sigma_Tc, sigma_V1, sigma_V2):
    # Compute X and Y
    ln_V = np.log(V1 / V2)
    X = R * (Th - Tc) * ln_V
    Y = cp * (Th - Tc) + R * Th * ln_V

    # Compute partial derivatives
    d_eta_dTh = (Y * R * ln_V - X * (cp + R * ln_V)) / Y**2
    d_eta_dTc = (Y * (-R * ln_V) - X * (-cp)) / Y**2
    d_eta_dV1 = (Y * (R * (Th - Tc) / V1) - X * (R * Th / V1)) / Y**2
    d_eta_dV2 = (Y * (-R * (Th - Tc) / V2) - X * (-R * Th / V2)) / Y**2

    print((d_eta_dTh * sigma_Th) ** 2 )
    print((d_eta_dTc * sigma_Tc) ** 2 )
    print((d_eta_dV1 * sigma_V1) ** 2 )
    print((d_eta_dV2 * sigma_V2) ** 2)

    # Compute overall uncertainty using propagation of errors
    sigma_eta = np.sqrt(
        (d_eta_dTh * sigma_Th) ** 2 +
        (d_eta_dTc * sigma_Tc) ** 2 +
        (d_eta_dV1 * sigma_V1) ** 2 +
        (d_eta_dV2 * sigma_V2) ** 2
    )

    return X/Y, sigma_eta

R = 287.101  # Specific gas constant for air (J/kg·K)
cp = 1005   # Specific heat capacity of air at constant pressure (J/kg·K)
Th = 273.15 + 46    # Hot temperature in K
Tc = 273.15 + 20    # Cold temperature in K
V1 = 61.38868395   # Initial volume in cm³
V2 = 52.91040895    # Final volume in cm³

# Uncertainties
sigma_Th = 0.2  # Uncertainty in Th (K)
sigma_Tc = 0.2  # Uncertainty in Tc (K)
sigma_V1 = 2.223733448 # Uncertainty in V1 (cm³)
sigma_V2 = 2.219690495 # Uncertainty in V2 (cm³)

# Compute uncertainty in eta
eta, sigma_eta = uncertainty_eta(R, cp, Th, Tc, V1, V2, sigma_Th, sigma_Tc, sigma_V1, sigma_V2)

print(f'Efficiency: {eta}')
print(f"Uncertainty in η: {sigma_eta}")
