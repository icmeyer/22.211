######################################################################
#               __       ____   __    __  __                         #
#              / /      /  _/  / /    \ \/ /                         #
#             / /       / /   / /      \  /                          #
#            / /___   _/ /   / /___    / /                           #
#           /_____/  /___/  /_____/   /_/                            #
#                                                                    #
######################################################################
#
#LILY is a 1D 1E Monte Carlo Code
#Author: Isaac Meyer
#Assumptions: -Isotropic scattering in COM
#             -Uniformly distributed isotropic unit source in slab 1
#             -Single energy group
#             -Vacuum boundary conditions
#Geometry: -Slab 1: 2cm
#          -Slab 2: 4cm
#Materials Propety: -Material 1: Total macro xs: 1 cm^-1
#                                Absor macro xs: 0.5 cm^-1
#                   -Material 2: Total macro xs: 1.5 cm^-1
#                                Absor macro xs: 1.2 cm^-1
# 
#Code structure: For simplicity I will be included all used functions
#                in this one file. I will not be using classes, but 
#                will emulate OOP in that all neutron attributes will
#                be variable names such as "neutron_e"  "neutron_xhat"
#                and so on
#
#Units: Distance: cm
#       Energy: eV 
#       Angle: radians
#
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import random_sample as rand
np.random.seed(42)

n_materials = 2
n_regions = 2
###################
### Create grid ###
###################
grid_width = 0.5
total_width = 6.0
n_cells = int(total_width/grid_width)
grid_boundaries = np.linspace(0,6,n_cells+1) 
#Assigning materials to grid
region_lengths = [2.0]
grid_mat = np.zeros([n_cells])
material_counter = 1
for i in range(n_cells):
    cell_start=grid_boundaries[i+1]
    if cell_start == region_lengths: #Make sure all region boundaries
                                     #align with cell boundaries
        material_counter += 1
    grid_mat[i] = material_counter
print(grid_mat)
    

