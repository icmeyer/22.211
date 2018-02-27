import numpy as np
from numpy.random import random_sample as rand
import matplotlib.pyplot as plt
np.random.seed(42)

def inv_polar(xi):
    mu=1/((16*xi**2-16*xi+5)**(1/2)-4*xi+2)**(1/3)-((16*xi**2-16*xi+5)**(1/2)-4*xi+2)**(1/3)
    return mu

trials = int(1e6)
mu_sum = 0
for i in range(trials):
    xi = rand()
    mu = inv_polar(xi)
    mu_sum += mu
print("Average:",mu_sum/trials)

