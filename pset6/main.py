######################################################################
#      _____         _____    ____     ______     ____     __        #
#     |    _ `''.  .'  __ `. \ .-.')  / _     \   \   \   /  /       #
#     | _ | ) _  \/   '  \  \/ `-' \ (`' )/`--'    \  _. /  '        #
#     |( ''_'  ) ||___|  /  | `-'`"`(_ o _).        _( )_ .'         #
#     | . (_) `. |   _.-`   | .---.  (_,_). '.  ___(_ o _)'          #
#     |(_    ._) '.'   _    | |   | .---.  \  :|   |(_,_)'           #
#     |  (_.\.' / |  _( )_  | |   | \    `-'  ||   `-'  /            #
#     |       .'  \ (_ o _) / |   |  \       /  \      /             #
#     '-----'`     '.(_,_).'  '---'   `-...-'    `-..-'              #
#                                                                    #
#     Diffusion Actually Isn't So Yucky                              #
######################################################################
from materials import mat_properties
import numpy as np
import matplotlib.pyplot as plt
from build_matrix import build_matrix
from power import power_iteration, normalize
from math import pi
"""
1D - 2 Group - Diffusion Equation Solver

This is the main script of the solver. It will have dependencies
on other files in order to run.
"""

#Initiate problem statements as a dictionary of dictionaries
#Each member of "problem_statement" is a separate problem
#Each member of a problem is a zone with list: 
#                           [Width [cm],Spacing,Material ID] 
problem_statement = {'1': {'1': [300,5.0,'1']},
                     '2': {'1': [25 ,5.0,'5'], '2': [250, 5.0, '2'],
                           '3': [25 ,5.0,'5']},
                     '3': {'1': [25 ,1.0,'5'], '2': [250, 1.0, '2'],
                           '3': [25 ,1.0,'5']},
                     '4': {'1': [25 ,1.0,'5'], '2': [15 , 1.0, '4'],
                           '3': [220,1.0,'3'], '4': [15 , 1.0, '4'],
                           '5': [25 ,1.0,'5']},
                     '5': {'1': [23 ,1.0,'7'], '2': [2  , 1.0, '6'],
                           '3': [15 ,1.0,'4'], '4': [220, 1.0, '3'],
                           '5': [15 ,1.0,'4'], '6': [2  , 1.0, '6'],
                           '7': [23 ,1.0,'7']},
                     'test': {'1': [30 ,5.0 ,'1']}
                     }


run = ['test']
run = ['1','2','3','4','5']
run = ['3']
for i in run:
    problem = i
    statement = problem_statement[problem]
    hmat, fmat , ncells= build_matrix(statement)
    # print(hmat)
    # print(fmat)
    plt.spy(hmat)
    plt.show()
    plt.spy(fmat)
    plt.show()
    phi_guess = np.ones([2*ncells])
    k_guess = 1
    phi,k,counter,b = power_iteration(hmat,fmat,phi_guess,k_guess,ncells)
    phi1 = phi[:ncells]
    phi2 = phi[ncells:]
    fissionsource = b[:ncells]+b[ncells:]
    peak=np.amax(fissionsource)/np.average(fissionsource)
    
    
    
    # print(hmat)
    # plt.spy(hmat)
    # plt.show()
    # print(fmat)
    # plt.spy(fmat)
    # plt.show()
    width = 300
    x = np.linspace(0,width,ncells)
    plt.plot(x,phi1)
    plt.plot(x,phi2)
    plt.legend(['Group 1','Group 2'])
    if problem == '1':
        # Analytic solution for problem 1
        width = 300
        x = np.linspace(0,width,ncells)
        bg = (pi/width)
        sigma_s12 =  mat_properties['1']['S12']
        sigma_a1  =  mat_properties['1']['A1'] 
        sigma_a2  =  mat_properties['1']['A2'] 
        sigma_r1  = sigma_a1+sigma_s12
        d1 = mat_properties['1']['D1']
        d2 = mat_properties['1']['D2']
        nf1 =  mat_properties['1']['NF1']
        nf2 =  mat_properties['1']['NF2']
        kanal = (nf1/(sigma_r1+d1*bg**2) +
                sigma_s12/(sigma_r1+d1*bg**2)*nf2/(sigma_a2+d2*bg**2))
        phi1anal = np.cos(pi*(x-width/2)/width)
        phi2anal = phi1anal*(sigma_s12/(d2*(pi/width)**2+sigma_a2))
        phitotanal = np.concatenate([phi1anal,phi2anal])
        #Normalizing
        phitotanal = normalize(phitotanal)
        phi1anal=phitotanal[ncells:]
        phi2anal=phitotanal[:ncells]
        print("Analytical keff is:", kanal)

        plt.plot(x,phi1anal,'.')
        plt.plot(x,phi2anal,'.')
        plt.legend(['Group 1 Numerical','Group 2 Numerical','Group 1 Analytical',
                    'Group 2 Analytical'])
    plt.title('k='+str(k)+' with '+str(counter)+' iterations')
    plt.show()
    plt.plot(x,fissionsource)
    plt.title('Fission Density, Peak='+str(peak))
    plt.show()

