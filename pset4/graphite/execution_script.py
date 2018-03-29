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
    all_tallies = pincellfunction(pitch, 0.711) # natural enrichment
    all_pitches = all_pitches.append(all_tallies, ignore_index=True)

# Criticality search for critical Dy2O3 concentation
enrichs = np.linspace(0.5,2,50)
all_enrichs = pd.DataFrame()
for enrich in enrichs:
    all_tallies = pincellfunction(8.58, enrich) #pitch is 15 cm
    all_enrichs = all_enrichs.append(all_tallies, ignore_index=True)

#Write results to csv
os.chdir('..')
all_pitches.to_csv('results_pitch.csv')
all_enrichs.to_csv('results_enrich.csv')
