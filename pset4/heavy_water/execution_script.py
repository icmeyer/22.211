from pincellfunction import pincellfunction
import openmc
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pandas as pd
import os, shutil

if not os.path.isdir('work'):
    os.mkdir('work')
os.chdir('work')
if not os.path.isdir('figures'):
    os.mkdir('figures')

#Range of pitches beginning at 2 cm
pitches = 2*np.logspace(0,1,50)

all_pitches = pd.DataFrame()
for pitch in pitches:
    all_tallies = pincellfunction(pitch, 0) #dy2o3_conc is 0
    all_pitches = all_pitches.append(all_tallies, ignore_index=True)

# Criticality search for critical Dy2O3 concentation
concs = np.linspace(0,0.02,50)
all_concs = pd.DataFrame()
for conc in concs:
    all_tallies = pincellfunction(10.86, conc) #pitch is 10.85 cm
    all_concs = all_concs.append(all_tallies, ignore_index=True)

#Write results to csv
os.chdir('..')
all_pitches.to_csv('results_pitch.csv')
all_concs.to_csv('results_conc.csv')
