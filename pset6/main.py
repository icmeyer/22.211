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
import matplotlib as mpl
from build_matrix import build_matrix
from power import power_iteration, normalize
from math import pi
import os
mpl.rcParams['figure.figsize'] = [12.0, 6.0]
"""
1D - 2 Group - Diffusion Equation Solver

This is the main script of the solver. It will have dependencies
on other files in order to run.
"""
runbaffles = True
gif = False

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


run = ['3']
run = ['test']
run = ['1','2','3','4','5']
for i in run:
    problem = i
    statement = problem_statement[problem]
    hmat, fmat , ncells= build_matrix(statement)
    # print(hmat)
    # print(fmat)
    plt.subplot(121)
    plt.spy(hmat)
    plt.title('H')
    plt.subplot(122)
    plt.spy(fmat)
    plt.title('F')
    plt.show()
    phi_guess = np.ones([2*ncells])
    k_guess = 1
    phi,k,counter,b = power_iteration(hmat,fmat,phi_guess,k_guess,ncells)
    phi1 = phi[:ncells] #1st half of vector is 1st group
    phi2 = phi[ncells:] #2nd half is 2nd group

    #Handling fission source
    fissionsource = b[:ncells] #Fissionsource only in first group
    fissionsource = b[np.nonzero(b)] #Average only over NF!=0
    fissionsource = fissionsource/np.average(fissionsource)
    #Peak is calculated as maximum when averaged over fuel
    peak=np.amax(fissionsource)
    print(fissionsource.shape)
    nfuelcells = fissionsource.shape[0]
    # Repad with zeros for plotting
    zeropad = int((ncells-nfuelcells)/2)
    fissionsource = np.concatenate([np.zeros(zeropad),fissionsource,
                                    np.zeros(zeropad)])
    print(np.average(fissionsource))
    
    width = 300
    x = np.linspace(0,width,ncells)
    plt.subplot(121)
    plt.plot(x,phi1)
    plt.plot(x,phi2)
    plt.xlabel('x position [cm]')
    plt.ylabel('Normalized Flux')
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
        extrapolated = width + 4*d1
        phi1anal = np.cos(pi*(x-width/2)/extrapolated)
        phi2anal = phi1anal*(sigma_s12/(d2*(pi/extrapolated)**2+sigma_a2))
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
    plt.subplot(122)
    plt.plot(x,fissionsource)
    plt.xlabel('x position [cm]')
    plt.ylabel('Fission Source (Average = 1)')
    plt.title('Fission Density, Peak='+str(peak))
    plt.show()


peaks = []
widths = []
keffs = []
if runbaffles == True:
    directory = './figs/'
    if not os.path.exists(directory):
            os.makedirs(directory)
    baf_change = np.linspace(0,23,47)
    for baf_change in baf_change:
        bafflewidth = str(2+baf_change)
        bafflewidthval = 2+baf_change
        widths.append(bafflewidthval)
        baffle = {'1': [23-baf_change ,0.5,'7'], '2': [2+baf_change , 0.5, '6'],
              '3': [15 ,0.5,'4'], '4': [220, 0.5, '3'],
              '5': [15 ,0.5,'4'], '6': [2+baf_change , 0.5, '6'],
              '7': [23-baf_change ,0.5,'7']}
        hmat, fmat , ncells= build_matrix(baffle)
        phi_guess = np.ones([2*ncells])
        k_guess = 1
        phi,k,counter,b = power_iteration(hmat,fmat,phi_guess,k_guess,ncells)
        keffs.append(k)
        phi1 = phi[:ncells]
        phi2 = phi[ncells:]
        #Handling fission source
        fissionsource = b[:ncells] #Fissionsource only in first group
        fissionsource = b[np.nonzero(b)] #Average only over NF!=0
        fissionsource = fissionsource/np.average(fissionsource)
        #Peak is calculated as maximum when averaged over fuel
        peak=np.amax(fissionsource)
        peaks.append(peak)
        nfuelcells = fissionsource.shape[0]
        # Repad with zeros for plotting
        zeropad = int((ncells-nfuelcells)/2)
        fissionsource = np.concatenate([np.zeros(zeropad),fissionsource,
                                        np.zeros(zeropad)])
        width = 300
        x = np.linspace(0,width,ncells)
        plt.subplot(121)
        plt.plot(x,phi1)
        plt.plot(x,phi2)
        plt.legend(['Group 1','Group 2'])
        plt.title('k='+str(k)+' with baffle width'+bafflewidth+' cm')
        plt.subplot(122)
        plt.plot(x,fissionsource)
        plt.title('Fission Density, Peak='+str(peak))
        if bafflewidthval<10:
            filename = 'baffle_0'+bafflewidth+'cm'
        else:
            filename = 'baffle_'+bafflewidth+'cm'
        plt.savefig(directory+filename+'.png')
        plt.close()
    if gif:
        os.system('cd '+directory+
                  '; convert -delay 20 -loop 0 *.png baffle.gif')
    plt.subplot(121)
    plt.plot(widths,peaks)
    plt.xlabel('Baffle Width (cm)')
    plt.ylabel('Peak of Fission Source (Average=1)')
    plt.subplot(122)
    plt.plot(widths,keffs)
    plt.xlabel('Baffle Width (cm)')
    plt.ylabel('k_eff')
    plt.show()

