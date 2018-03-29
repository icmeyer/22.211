import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('results_pitch_square.csv')
plt.errorbar(df['mf_ratio'],df['kinf mean'],yerr=df['kinf sd'],
        fmt='.',markersize=12,capsize=3)
df = pd.read_csv('results_pitch_triangle.csv')
plt.errorbar(df['mf_ratio'],df['kinf mean'],yerr=df['kinf sd'],
        fmt='.',markersize=12,capsize=3)
plt.legend(['Square Lattice','Triangular Lattice'])
plt.xlabel('Moderator to Fuel Ratio')
plt.ylabel('$k_{inf}$')
plt.show()
