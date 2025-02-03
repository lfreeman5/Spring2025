import numpy as np
from scipy.optimize import fsolve
from scipy.interpolate import interp1d

T_KP_water_gas = np.array([500,1000,1500,2000,2500,3000,3500])
KP_water_gas = np.array([138.3,1.443,0.3887,0.22,0.1635,0.1378,0.1241])

T = np.array([200, 298, 300, 400, 500, 600, 700, 800, 900, 1000, 
              1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 
              2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 
              3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000])

O2_delta_h = np.array([-2836, 0, 54, 3031, 6097, 9254, 12503, 15838, 19250, 22721, 
                    26232, 29775, 33350, 36355, 40590, 44253, 47943, 51660, 55342, 
                    59169, 62959, 66773, 70609, 74467, 78346, 82245, 86164, 90003, 
                    94060, 98036, 102029, 106040, 110068, 114112, 118173, 122249, 
                    126341, 130448, 134570, 138705])

N2_delta_h = np.array([-2841, 0, 54, 2973, 5920, 8905, 11942, 15046, 18222, 21468, 
                       24770, 28118, 31510, 34939, 38404, 41999, 45423, 48971, 52541, 
                       56130, 59738, 63360, 66997, 70645, 74305, 77974, 81652, 85338, 
                       89031, 92730, 96436, 100148, 103865, 107587, 111315, 115048, 
                       118786, 122528, 126276, 130028])

H2O_delta_h = np.array([-3227, 0, 62, 3458, 6947, 10528, 14209, 18005, 21930, 25993, 
                        30191, 34181, 38963, 43520, 48181, 52939, 57786, 62717, 67725, 72805, 
                        77952, 83160, 88426, 93744, 99112, 104524, 109979, 115472, 121001, 126563, 
                        132156, 137777, 143426, 149099, 154795, 160514, 166252, 172011, 177787, 183582])

CO2_delta_h = np.array([-3423, 0, 69, 4003, 8301, 12899, 17749, 22810, 28047, 33425, 
                        38911, 44488, 50149, 55882, 61681, 67538, 73446, 79399, 85392, 91420, 
                        97477, 103562, 109670, 115798, 121944, 128107, 134284, 140474, 146677, 152891, 
                        159116, 165351, 171597, 177853, 184120, 190397, 196685, 202983, 209293, 215613])

CO_delta_h = np.array([-2835, 0, 54, 2979, 5943, 8955, 12029, 15176, 18401, 21697, 
                       25046, 28440, 31874, 35345, 38847, 42379, 45937, 49517, 53118, 56737, 
                       60371, 64020, 67682, 71354, 75036, 78727, 82426, 86132, 89844, 93562, 
                       97287, 101016, 104751, 108490, 112235, 115985, 119739, 123499, 127263, 131032])

O2_dh = interp1d(T, O2_delta_h, kind='cubic', fill_value="extrapolate")
N2_dh = interp1d(T, N2_delta_h, kind='cubic', fill_value="extrapolate")
H2O_dh = interp1d(T, H2O_delta_h, kind='cubic', fill_value="extrapolate")
CO2_dh = interp1d(T, CO2_delta_h, kind='cubic', fill_value="extrapolate")
CO_dh = interp1d(T, CO_delta_h, kind='cubic', fill_value="extrapolate")

def solve_wgs_equilibrium(T):
    interpKP = interp1d(T_KP_water_gas, KP_water_gas, kind='linear', fill_value="extrapolate")
    Kp = interpKP(T)
    def equations(b):
        c = 1 - b
        d = 14/6 - b
        e = 2 - d
        return Kp - (b * e) / (c * d)

    b_initial = 0.5  
    b_solution = fsolve(equations, b_initial)[0]
    c_solution = 1 - b_solution
    d_solution = 14/6 - b_solution
    e_solution = 2 - d_solution

    assert b_solution>1e-15
    return b_solution, c_solution, d_solution, e_solution

def delta_H(T):
    # Calculate the difference in enthalpy between reactants and products
    a=10./6.
    b,c,d,e = solve_wgs_equilibrium(T)

    H_reac = -74631
    H_prod = b*(-393546 + CO2_dh(T)) + c*(-110541 + CO_dh(T)) + d*(-241845 + H2O_dh(T)) + e*O2_dh(T) + (a*3.76)*N2_dh(T)
    print(f'Temperature guess: {T}K  Difference in enthalpy (Products-Reactants): {H_prod-H_reac}J/mol')
    return H_prod-H_reac

if __name__ == "__main__":
    T_guess = 2000.0
    
    T_solution = fsolve(delta_H, T_guess, full_output=True, xtol=1e-8)

    # Print the results
    print(f"Final temperature at equilibrium: {T_solution[0][0]:.2f} K")
    print(f'b,c,d,e values at equilibrium: {solve_wgs_equilibrium(T_solution[0][0])}')
    b,c,d,e = solve_wgs_equilibrium(T_solution[0][0])
    N_tot = 10./6.*3.76 +b+c+d+e
    print(N_tot)
    print(f'Xco2:{b/N_tot}')
    print(f'Xco:{c/N_tot}')
    print(f'XH20:{d/N_tot}')
    print(f'XH2:{e/N_tot}')
    print(f'XN2:{10./6.*3.76/N_tot}')
    