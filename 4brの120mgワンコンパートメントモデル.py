import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

def calculate_concentration_curves(params):
    t = np.arange(0, 72, 0.1)
    ka = (params['F'] * params['DOSE'] / params['Vd']) / (1 + 10 ** ((params['pH'] - 6.56) * np.log10(np.exp(1)) / np.log10(10)) * np.exp(-1 * (params['F'] * params['DOSE'] / params['Vd']) * params['onset'] / params['BW']))
    ke = np.log(2) / (params['duration'] + params['aftereffect'])
    Cp = (ka / (ka - ke)) * (np.exp(-1 * ke * t) - np.exp(-1 * ka * t)) * params['DOSE'] / params['Vd']
    return t, Cp

def plot_concentration_curve(params):
    t, Cp = calculate_concentration_curves(params)
    AUC = simps(Cp, t)
    
    plt.plot(t, Cp, label='Blood Concentration', color='blue')
    plt.ylabel('Blood Concentration (mg/L)')
    plt.ylim(0, 120)
    plt.xlabel('Time (h)')
    plt.xlim(0, 72)
    plt.title(f"Blood Concentration Curve with Dose: {params['DOSE']} mg")
    plt.legend(loc='upper left')
    plt.text(40, 100, f'AUC(h*mg/L): {AUC:.2f}', fontsize=12, bbox=dict(facecolor='red', alpha=0.5))
    plt.tight_layout()
    plt.show()

params = {
    'DOSE': 120,
    'BW': 60,
    'onset': 1/60,
    'duration': 12,
    'aftereffect': 28,
    'pH': 7.4,
    'F': 0.95,
    'Vd': 3.35
}

plot_concentration_curve(params)
