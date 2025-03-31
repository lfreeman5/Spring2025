import numpy as np
import math

# Given values
h_w = 2257  # kJ/kg (constant)

# Lists of values
HR1 = [0.006, 0.006, 0.006, 0.0065, 0.0065]
HR2 = [0.0064, 0.0062, 0.0069, 0.0213, 0.0201]
h1 = [37.63159, 37.63159, 37.63159, 39.72609, 39.72609]
h2 = [38.34, 38.33, 56.04, 96.36, 77.41]
h3 = [36.95, 38.3, 53.46, 93.77, 77.36]
h4 = [35.58, 48.4, 60.08, 98.35, 86]

# Flow rates calculated previously (assuming these are the values from the previous calculation)
# You'll need to replace these with the actual values from your previous calculation
m_dot_2 = 6*np.array([0.028026,0.026802,0.025844,0.025195,0.025855])  # Example values, replace with your calculated values
sigma_m_dot_2 = [0.000441,0.000415,0.000394,0.000395,0.000411]  # Example values, replace with your calculated uncertainties

# Uncertainties
sigma_h = 3.0  # For all h values
sigma_HR = 0.0004  # For all HR values

# Function to calculate Qa, Qb, Qc
def calculate_Q_values(m_dot, hr1, hr2, h1_val, h2_val, h3_val, h4_val, h_w):
    # Calculate term for both Qa and Qb
    term1 = (1 + hr1) / (1 + hr2)
    term2 = (hr2 - hr1) / (1 + hr2)
    
    # Calculate Qb
    Qb = m_dot * (h2_val - term1 * h1_val - term2 * h_w)
    
    # Calculate Qc
    Qc = m_dot * (h4_val - h3_val)
    
    # Calculate Qa
    Qa = Qb + Qc
    
    return Qa, Qb, Qc

# Function to calculate uncertainty of Qb
def calculate_uncertainty_Qb(m_dot, sigma_m_dot, hr1, hr2, h1_val, h2_val, h_w, sigma_hr, sigma_h):
    # Common terms
    term1 = (1 + hr1) / (1 + hr2)
    term2 = (hr2 - hr1) / (1 + hr2)
    
    # Term for m_dot uncertainty
    term_m_dot = (h2_val - term1 * h1_val - term2 * h_w) * sigma_m_dot
    
    # Term for HR1 uncertainty
    term_hr1 = m_dot * (h1_val / (1 + hr2)**2 - h_w / ((1 + hr2))) * sigma_hr
    
    # Term for HR2 uncertainty
    term_hr2 = m_dot * (-((1 + hr1) * h1_val) / (1 + hr2)**2 + h_w / (1 + hr2)**2) * sigma_hr
    
    # Term for h1 uncertainty
    term_h1 = m_dot * (-(1 + hr1) / (1 + hr2)) * sigma_h
    
    # Term for h2 uncertainty
    term_h2 = m_dot * sigma_h
    
    # Calculate final uncertainty
    sigma_Qb = math.sqrt(term_m_dot**2 + term_hr1**2 + term_hr2**2 + term_h1**2 + term_h2**2)
    
    return sigma_Qb

# Function to calculate uncertainty of Qc
def calculate_uncertainty_Qc(m_dot, sigma_m_dot, h3_val, h4_val, sigma_h):
    # Term for m_dot uncertainty
    term_m_dot = (h4_val - h3_val) * sigma_m_dot
    
    # Term for h4 uncertainty
    term_h4 = m_dot * sigma_h
    
    # Term for h3 uncertainty
    term_h3 = m_dot * sigma_h
    
    # Calculate final uncertainty
    sigma_Qc = math.sqrt(term_m_dot**2 + term_h4**2 + term_h3**2)
    
    return sigma_Qc

# Calculate and store results
Qa_values = []
Qb_values = []
Qc_values = []
sigma_Qb_values = []
sigma_Qc_values = []

for i in range(len(m_dot_2)):
    Qa, Qb, Qc = calculate_Q_values(
        m_dot_2[i], HR1[i], HR2[i], h1[i], h2[i], h3[i], h4[i], h_w
    )
    
    sigma_Qb = calculate_uncertainty_Qb(
        m_dot_2[i], sigma_m_dot_2[i], HR1[i], HR2[i], h1[i], h2[i], h_w, sigma_HR, sigma_h
    )
    
    sigma_Qc = calculate_uncertainty_Qc(
        m_dot_2[i], sigma_m_dot_2[i], h3[i], h4[i], sigma_h
    )
    
    Qa_values.append(Qa)
    Qb_values.append(Qb)
    Qc_values.append(Qc)
    sigma_Qb_values.append(sigma_Qb)
    sigma_Qc_values.append(sigma_Qc)

# Display results
print("Results for heat transfer calculations:")
print("------------------------------------------------------")
print("Trial  Qa (kJ/s)  Qb (kJ/s)  Qc (kJ/s)  Uncertainty Qb (kJ/s)  Uncertainty Qc (kJ/s)")

for i in range(len(m_dot_2)):
    print(f"{i+1}  {Qa_values[i]:.4f}  {Qb_values[i]:.4f}  {Qc_values[i]:.4f}  {sigma_Qb_values[i]:.4f}  {sigma_Qc_values[i]:.4f}")

# Calculate averages
avg_Qa = sum(Qa_values) / len(Qa_values)
avg_Qb = sum(Qb_values) / len(Qb_values)
avg_Qc = sum(Qc_values) / len(Qc_values)

# For uncertainty propagation in the average, we use the root sum of squares method
avg_sigma_Qb = math.sqrt(sum(u**2 for u in sigma_Qb_values)) / len(sigma_Qb_values)
avg_sigma_Qc = math.sqrt(sum(u**2 for u in sigma_Qc_values)) / len(sigma_Qc_values)

print("------------------------------------------------------")
print(f"Average Qa: {avg_Qa:.4f} kJ/s")
print(f"Average Qb: {avg_Qb:.4f} ± {avg_sigma_Qb:.4f} kJ/s")
print(f"Average Qc: {avg_Qc:.4f} ± {avg_sigma_Qc:.4f} kJ/s")

# Print copy-paste friendly output
print("\n# Copy-paste friendly output:")
print("Qa_values = [", end="")
for i, val in enumerate(Qa_values):
    if i < len(Qa_values) - 1:
        print(f"{val:.4f}, ", end="")
    else:
        print(f"{val:.4f}]")

print("Qb_values = [", end="")
for i, val in enumerate(Qb_values):
    if i < len(Qb_values) - 1:
        print(f"{val:.4f}, ", end="")
    else:
        print(f"{val:.4f}]")

print("Qc_values = [", end="")
for i, val in enumerate(Qc_values):
    if i < len(Qc_values) - 1:
        print(f"{val:.4f}, ", end="")
    else:
        print(f"{val:.4f}]")

print("sigma_Qb_values = [", end="")
for i, val in enumerate(sigma_Qb_values):
    if i < len(sigma_Qb_values) - 1:
        print(f"{val:.4f}, ", end="")
    else:
        print(f"{val:.4f}]")

print("sigma_Qc_values = [", end="")
for i, val in enumerate(sigma_Qc_values):
    if i < len(sigma_Qc_values) - 1:
        print(f"{val:.4f}, ", end="")
    else:
        print(f"{val:.4f}]")

print(f"avg_Qa = {avg_Qa:.4f}")
print(f"avg_Qb = {avg_Qb:.4f}")
print(f"avg_Qc = {avg_Qc:.4f}")
print(f"avg_sigma_Qb = {avg_sigma_Qb:.4f}")
print(f"avg_sigma_Qc = {avg_sigma_Qc:.4f}")