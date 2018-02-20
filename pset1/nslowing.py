import numpy as np
import matplotlib.pyplot as plt
from math import log
def nslowing(bins):
    neutrons = 10000
    e_bins = np.linspace(0,100,bins)
    leth_bins = np.logspace(0,2,bins)
    e_freq = np.zeros(bins-1)
    leth_freq = e_freq
    collision_dist = {1: e_freq,2: e_freq,3: e_freq}
    
    for i in range(neutrons):
        energy = 1000 
        collisions = 0
        while energy > 1:
            e_freq_step = np.histogram(energy,e_bins)[0]
            e_freq = e_freq+e_freq_step/20 #Divide by background xs for flux
            leth_freq_step = np.histogram(energy,leth_bins)[0]
            leth_freq = leth_freq+leth_freq_step/20
            if collisions in [1,2,3]:
                collision_dist[collisions] = collision_dist[collisions]+e_freq_step
            energy = np.random.uniform(0,1)*energy
            collisions = collisions+1
    
    e_freq = np.concatenate([[0],e_freq])/neutrons
    leth_freq = np.concatenate([[0],leth_freq])/neutrons
    for i in range(3):
        collision_dist[i+1] = np.concatenate([[0],collision_dist[i+1]])/neutrons
    return e_bins, leth_bins, e_freq, leth_freq, collision_dist
