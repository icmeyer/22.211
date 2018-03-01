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
#Materials Propety: -Material 0: Total macro xs: 1 cm^-1
#                                Absor macro xs: 0.5 cm^-1
#                   -Material 1: Total macro xs: 1.5 cm^-1
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
import math
from math import pi,sin,cos
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
grid_mats = np.zeros([n_cells])
material_counter = 0
for i in range(n_cells):
    cell_start=grid_boundaries[i+1]
    grid_mats[i] = material_counter
    if cell_start == region_lengths: #Make sure all region boundaries
                                     #align with cell boundaries
        material_counter += 1
print(grid_boundaries)
print(grid_mats)

#################
# Material Data #
#################
#Material 0: Total macro xs: 1 cm^-1
#            Absor macro xs: 0.5 cm^-1
#Material 1: Total macro xs: 1.5 cm^-1
#            Absor macro xs: 1.2 cm^-1
def get_grid_info(grid_boundaries,grid_mats,neutron_x):
    # Constants for problem
    mats_xs_total = [1,1.5]
    mats_xs_abs = [0.5,1.2]
    ################################################3
    grid_location_bool = np.histogram(neutron_x,grid_boundaries)[0]
    try:
        cell = int(np.where(grid_location_bool)[0]) #Returns the cell number
        material = int(grid_mats[cell])
        xs_total = mats_xs_total[material]
        xs_abs = mats_xs_abs[material]
    except: #leaked
        cell = []; material = []; xs_total = []; xs_abs=[];
    return xs_total,xs_abs,cell,material

    
#############
# Functions #
#############
def sample_d(macro_xs):
    x = -math.log(rand())/macro_xs
    return x

def rand_unit_vector():
    vect = [2*rand()-1,2*rand()-1,2*rand()-1]
    vect = vect/np.linalg.norm(vect)
    return vect

def rand_mu_l():
    A = 1
    mu_cm = 2*rand()-1
    mu_l = (1+A*mu_cm)/(A**2+2*A*mu_cm+1)**(1/2)
    return mu_l
    
def rhat_prime(rhat_naut):
    mu_l = rand_mu_l()
    phi = 2*pi*rand()
    u = rhat_naut[0]; v=rhat_naut[1]; w=rhat_naut[2];
    uprime = mu_l*u+(1-mu_l**2)**(1/2)*(u*w*cos(phi)-v*sin(phi))/(1-w**2)**(1/2)
    vprime = mu_l*v+(1-mu_l**2)**(1/2)*(v*w*cos(phi)+u*sin(phi))/(1-w**2)**(1/2)
    wprime = mu_l*w-(1-mu_l**2)**(1/2)*(1-w**2)**(1/2)*cos(phi)
    return uprime,vprime,wprime

#########################
### Run Some Neutrons ###
#########################
# Initialize neutron
nps = 1000
for i in range(nps):
    print("~~~~~~~~~~~~~~ Baby Neutron:",i,"~~~~~~~~~~~~~~~")
    ##### Initialize our baby neutron #####
    neutron_x = rand()*total_width
    neutron_rhat = rand_unit_vector()
    neutron_v = 0 #Unused in this problem
    while True:
        print("Position",neutron_x,"Unit Vector",neutron_rhat,np.linalg.norm(neutron_rhat))
        xs_total,xs_abs,cell,material = get_grid_info(grid_boundaries,grid_mats,neutron_x)
        ##### Move that baby #####
        neutron_x += neutron_rhat[0]*sample_d(xs_total)
        print(neutron_x)
        xs_total,xs_abs,new_cell,new_material = get_grid_info(grid_boundaries,grid_mats,neutron_x)

        ##### Handling leakage and change of cell #####
        if not new_cell:
            #This shit totally leaked
            print("Shit totally leaked")
            break
        elif new_cell != cell:
            print("changed cell")
            while True:
                print("loop")
                if neutron_rhat[0] < 0: #Came from the right
                    cell = cell-1
                    neutron_x = grid_boundaries[cell+1]+neutron_rhat[0]*sample_d(xs_total)
                    xs_total,xs_abs,new_cell,new_material = get_grid_info(grid_boundaries,grid_mats,neutron_x)
                    if new_cell == cell:
                        break
                    elif not new_cell:
                        print("Shit totally leaked")
                        break
                else: #Came from the left
                    cell = cell+1
                    neutron_x = grid_boundaries[cell]+neutron_rhat[0]*sample_d(xs_total)
                    xs_total,xs_abs,new_cell,new_material = get_grid_info(grid_boundaries,grid_mats,neutron_x)
                    if new_cell == cell:
                        break
                    elif not new_cell:
                        print("Shit totally leaked")
                        break
        print("Position",neutron_x,"Unit Vector",neutron_rhat,np.linalg.norm(neutron_rhat))
        if not new_cell:
            #This shit totally leaked
            print("Shit totally leaked")
            break

        ##### Actual collision #####
        #Sample isotope (not present because homegeneous in this case
        #Sample collision type
        if rand()<xs_abs/xs_total:
            print("Baby got eaten")
            break
        else:
            print("Baby got bounced!")
            neutron_rhat = rhat_prime(neutron_rhat)

