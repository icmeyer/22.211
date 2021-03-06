import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('results_pitch.csv')
plt.errorbar(df['pitch'],df['kinf mean'],yerr=df['kinf sd'],
        fmt='.',markersize=12,capsize=3)
plt.xlabel('Pitch [cm]')
plt.ylabel('$k_{inf}$')
plt.title('Natural Enrichment, Heavy Water Moderator')
plt.show()

plt.errorbar(df['pitch'],df['res_esc mean'],yerr=df['res_esc sd'],
        fmt='.',markersize=12,capsize=3)
plt.errorbar(df['pitch'],df['fast_fiss mean'],yerr=df['fast_fiss sd'],
        fmt='.',markersize=12,capsize=3)
plt.errorbar(df['pitch'],df['therm_util mean'],yerr=df['therm_util sd'],
        fmt='.',markersize=12,capsize=3)
plt.errorbar(df['pitch'],df['eta mean'],yerr=df['eta sd'],
        fmt='.',markersize=12,capsize=3)
plt.legend(['Resonance Escape, $p$','Fast Fission, $\epsilon$','Thermal Utilization $f$','Reproduction Factor, $\eta$'])
plt.xlabel('Pitch [cm]')
plt.title('Four Factor Values: Natural Enrichment, Heavy Water Moderator')
plt.show()

df = pd.read_csv('results_conc.csv')
plt.errorbar(df['dy2o3 conc'],df['kinf mean'],yerr=df['kinf sd'],
        fmt='.',markersize=12,capsize=3)
plt.xlabel('$Dy_2 O_3$ conc. [parts per $UO_2$]')
plt.ylabel('$k_{inf}$')
plt.title('Natural Enrichment, Heavy Water Moderator, Pitch = 11.38 cm')
plt.show()
