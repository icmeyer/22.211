import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('results_pitch.csv')
plt.errorbar(df['pitch'],df['kinf mean'],yerr=df['kinf sd'],
        fmt='.',markersize=12,capsize=3)
plt.xlabel('Pitch [cm]')
plt.ylabel('$k_{inf}$')
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
plt.title('Four Factor Formula Values')
plt.show()

plt.errorbar(df['pitch'],df['res_esc mean'],yerr=df['res_esc sd'],
        fmt='.',markersize=12,capsize=3)
plt.errorbar(df['pitch'][4:],df['fast_fiss mean'][4:],yerr=df['fast_fiss sd'][4:],
        fmt='.',markersize=12,capsize=3)
plt.errorbar(df['pitch'],df['therm_util mean'],yerr=df['therm_util sd'],
        fmt='.',markersize=12,capsize=3)
plt.errorbar(df['pitch'],df['eta mean'],yerr=df['eta sd'],
        fmt='.',markersize=12,capsize=3)
plt.legend(['Resonance Escape, $p$','Fast Fission, $\epsilon$','Thermal Utilization $f$','Reproduction Factor, $\eta$'])
plt.xlabel('Pitch [cm]')
plt.title('Four Factor Formula Values')
plt.show()

df = pd.read_csv('results_conc.csv')
plt.errorbar(df['enrichment'],df['kinf mean'],yerr=df['kinf sd'],
        fmt='.',markersize=12,capsize=3)
plt.xlabel('Uranium Enrichment [weight percent]')
plt.ylabel('$k_{inf}$')
plt.show()
