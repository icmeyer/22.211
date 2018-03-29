from pincellfunction import pincellfunction
from math import pi
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
    all_tallies = pincellfunction(pitch, 4.0) #enrichment is 4.0 percent
    all_pitches = all_pitches.append(all_tallies, ignore_index=True)

#k as a function of volume ratio
pitches = all_pitches['pitch']
radius = 0.5
ratio = ((3**(1/2)/4)*pitches**2-pi*radius**2)/(pi*radius**2)
all_pitches = all_pitches.assign(mf_ratio = ratio)

#Write results to csv
os.chdir('..')
all_pitches.to_csv('results_pitch.csv')
