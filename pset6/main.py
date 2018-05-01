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
import materials
import numpy as np
import matplotlib.pyplot as plt
from build_matrix import build_matrix
from power import power_iteration
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
                     'test': {'1': [3,1.0,'1']}
                     }

problem = problem_statement['5']
hmat, fmat , ncells= build_matrix(problem)
phi_guess = np.ones([2*ncells])/(2*ncells)
k_guess = 1
phi, k, counter = power_iteration(hmat,fmat,phi_guess,k_guess,ncells)
phi1 = phi[:ncells]
phi2 = phi[ncells:]

# print(hmat)
# plt.spy(hmat)
# plt.show()
# print(fmat)
# plt.spy(fmat)
# plt.show()
x = np.linspace(0,300,ncells)
plt.plot(x,phi1)
plt.plot(x,phi2)
plt.title('k='+str(k))
plt.show()

