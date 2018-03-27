from pincellfunction import pincell
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
pitches = 2*np.logspace(0,2,5)

all_pitches = pd.DataFrame()
for pitch in pitches:
    all_tallies = pincell(pitch)
    print(all_tallies)
    all_pitches.append(all_tallies, ignore_index=True)

#Write results to csv
os.chdir('..')
all_pitches.to_csv('results.csv')



