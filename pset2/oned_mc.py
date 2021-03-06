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
grid_width = .25
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

###################
# Create tallies  #
###################
collision_based_s1 = np.zeros([n_cells])
collision_based_s2 = np.zeros([n_cells])
path_length_s1 = np.zeros([n_cells])
path_length_s2 = np.zeros([n_cells])

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
        cell = 'void'; material = []; xs_total = []; xs_abs=[];
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
nps = 10000
for i in range(nps):
    print("~~~~~~~~~~~~~~ Baby Neutron:",i,"~~~~~~~~~~~~~~~")
    ##### Initialize counters #####
    collision_counter = np.zeros([n_cells]) 
    path_length_counter = np.zeros([n_cells])
    ##### Initialize our baby neutron #####
    neutron_x = rand()*region_lengths[0]
    neutron_rhat = rand_unit_vector()
    neutron_v = 0 #Unused in this problem
    while True:
        xs_total,xs_abs,cell,material = get_grid_info(grid_boundaries,grid_mats,neutron_x)
        print("Before Distance Sample: Position",neutron_x,"Unit Vector",neutron_rhat,"Cell",cell)
        ##### Move that baby #####
        path_length = neutron_rhat[0]*sample_d(xs_total)
        new_neutron_x = neutron_x + path_length
        new_xs_total,new_xs_abs,new_cell,new_material = get_grid_info(grid_boundaries,grid_mats,new_neutron_x)

        ##### Handling leakage and change of cell #####
        if new_cell != cell:
            print("Sampled x in new cell. New x:",new_neutron_x)
            while True: #This loop will run until the neutron settles in a cell
                if neutron_rhat[0] < 0: #Came from the right
                    path_length_counter[cell] += neutron_x-grid_boundaries[cell]
                    cell = cell-1
                    if cell == -1:
                        new_cell = 'void'
                        print("Leaked to left")
                        break
                    else:
                        neutron_x = grid_boundaries[cell+1]
                        new_xs_total = get_grid_info(grid_boundaries,grid_mats,neutron_x-0.0001)[0]
                        path_length = neutron_rhat[0]*sample_d(new_xs_total)
                        new_neutron_x = neutron_x + path_length
                        new_xs_total,new_xs_abs,new_cell,new_material = get_grid_info(grid_boundaries,grid_mats,new_neutron_x)
                    if new_cell == cell: #Break loop if neutron stayed in same cell
                        break
                else: #Came from the left
                    path_length_counter[cell] += grid_boundaries[cell+1]-neutron_x
                    cell = cell+1
                    if cell == n_cells:
                        new_cell = 'void'
                        print("Leaked to right")
                        break
                    else:
                        neutron_x = grid_boundaries[cell]
                        new_xs_total = get_grid_info(grid_boundaries,grid_mats,neutron_x+0.0001)[0]
                        path_length = neutron_rhat[0]*sample_d(new_xs_total)
                        new_neutron_x = neutron_x + path_length
                        new_xs_total,new_xs_abs,new_cell,new_material = get_grid_info(grid_boundaries,grid_mats,new_neutron_x)
                    if new_cell == cell: #Break loop if neutron stayed in same cell
                        break

        if new_cell=='void': #Remove leaked neutrons
            print("Shit totally leaked: path_length_counter:",path_length_counter)
            collision_based_s1 += collision_counter
            collision_based_s2 += collision_counter**2
            path_length_s1 += path_length_counter
            path_length_s2 += path_length_counter**2
            break #Exits major loop
        if new_cell == cell: #Path length counter for movement in same cell
            neutron_x = new_neutron_x
            path_length_counter[cell] += abs(path_length)

        ##### Actual collision #####
        collision_counter[new_cell] += 1/new_xs_total
        print("Before collision: Position",new_neutron_x,"Unit Vector",neutron_rhat,"Cell",new_cell)
        print(collision_counter)
        #Sample isotope (not present because homegeneous in this case
        #Sample collision type
        if rand()<new_xs_abs/new_xs_total:
            print("Baby got eaten: path_length_counter:",path_length_counter)
            collision_based_s1 += collision_counter
            collision_based_s2 += collision_counter**2
            path_length_s1 += path_length_counter
            path_length_s2 += path_length_counter**2
            break #Exits major loop
        else:
            print("Baby got bounced!")
            neutron_rhat = rhat_prime(neutron_rhat)
        
collision_based_mean = collision_based_s1/sum(collision_based_s1)
collision_based_sd = ((1/(nps-1))*(collision_based_s2/nps-collision_based_mean**2))**(1/2)
path_length_mean = path_length_s1/sum(path_length_s1)
path_length_sd = ((1/(nps-1))*(path_length_s2/nps-path_length_mean**2))**(1/2)

print('Collision')
print(np.vstack([collision_based_mean,collision_based_sd]))
print('Path Length')
print(np.vstack([path_length_mean,path_length_sd]))
print('Divided')
print(collision_based_mean/path_length_mean)
print('Divided after normalization')
print(collision_based_mean/sum(collision_based_mean)/path_length_mean/sum(path_length_mean))
print(grid_boundaries)

#### Plotting Junk ####
collision_plot = np.concatenate([[0],collision_based_mean])
pl_plot = np.concatenate([[0],path_length_mean])
plt.step(grid_boundaries,collision_plot,color='C0')
plt.errorbar(grid_boundaries[:-1]+0.5*grid_width,collision_based_mean,yerr=collision_based_sd,fmt='o',markersize=1,capsize=3,color='C0')
plt.step(grid_boundaries,pl_plot,color='C1')
plt.errorbar(grid_boundaries[:-1]+0.5*grid_width,path_length_mean,yerr=path_length_sd,fmt='o',markersize=1,capsize=3,color='C1')
plt.legend(['Collison Based','Path Length Based'])
plt.xlabel('1d Cells (cm)')
plt.ylabel('Relative Flux (dimensionless for 1D Problem)')
plt.show()


