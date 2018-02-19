import numpy as np
import matplotlib.pyplot as plt
from math import log

neutrons = 10000
bins = 100
e_bins = np.linspace(1,100,bins)
leth_bins = np.logspace(0,2,bins)
e_freq = np.zeros(bins-1)
leth_freq = e_freq
collision_dist = {1: e_freq,2: e_freq,3: e_freq}

for i in range(neutrons):
    energy = 1000 
    collisions = 0
    while energy > 1:
        e_freq_step = np.histogram(energy,e_bins)[0]
        e_freq += e_freq_step
        leth_freq_step = np.histogram(energy,leth_bins)[0]
        leth_freq += leth_freq_step
        if collisions in [1,2,3]:
            collision_dist[collisions] = collision_dist[collisions]+e_freq_step
        energy = np.random.uniform(0,1)*energy
        collisions = collisions+1

e_freq = np.concatenate([[0],e_freq])/neutrons
ax = plt.subplot(3,1,1)
plt.step(e_bins,e_freq)
plt.ylabel('Normalized Flux')
plt.title('Equal Energy Bins n='+str(neutrons))
leth_freq = np.concatenate([[0],leth_freq])/neutrons
ax = plt.subplot(3,1,2)
ax.set_xscale("log", nonposx='clip')
plt.step(leth_bins,leth_freq)
plt.ylabel('Normalized Flux')
plt.title('Equal Lethargy Bins n='+str(neutrons))
ax=plt.subplot(3,1,3)
ax.set_xscale("log", nonposx='clip')
for i in range(3):
    collision_dist[i+1] = np.concatenate([[0],collision_dist[i+1]])/neutrons
    plt.step(leth_bins,collision_dist[i+1])
plt.legend(['1 Collison','2 Collisions','3 Collisions'])
plt.title('Energy Distribution after collisions (Equal Lethargy)')
plt.xlabel('Energy [eV]')
plt.ylabel('Normalized Flux')

plt.show()


