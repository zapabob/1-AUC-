import numpy as np
import matplotlib.pyplot as plt

# Parameters
dose = 180 # mg
half_life = 2  #hours
duration = 12 #hours
dt = 0.01  # hours
weight = 90  # kg
pH = 7.4
Vd = 20 

# Calculate clearance rate
clearance_rate = 0.693 / half_life

# Simulate drug concentration over time
times = np.arange(0, duration + dt, dt)
concentrations = np.zeros_like(times)
concentrations[0] = dose / Vd
for i in range(1, len(times)):
    concentrations[i] = concentrations[i-1] * np.exp(-clearance_rate * dt)

# Calculate the AUC
auc = np.trapz(concentrations, dx=dt)

# Plot results
plt.plot(times, concentrations, label='Drug Concentration')
plt.axhline(y=concentrations[-1], color='r', linestyle='--',label='Steady state')
plt.fill_between(times, concentrations, 0, alpha=0.1)
plt.xlabel('Time [h]')
plt.ylabel('Drug concentration [mg/mL]')
plt.title('AUC Curve')
plt.text(0, 1.1*concentrations[0], f'AUC: {auc:.2f} (mg*h/mL)')
plt.legend()
plt.show()
plt.savefig('auc_curve.png', dpi=300)
