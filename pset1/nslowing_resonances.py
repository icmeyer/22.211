import numpy as np
import matplotlib.pyplot as plt
from math import log
from res_construction import xs_from_res


def zerod_mc(mf_ratio,temp, neutrons, logemin, logemax, nbins, pick_res=[False,0]):
    #Resonances used in this problem
    res_E = [6.673491e+0, 2.087152e+1, 3.668212e+1]
    J = [5.000000e-1, 5.000000e-1, 5.000000e-1]
    gn = [1.475792e-3, 1.009376e-2, 3.354568e-2]
    gg = [2.300000e-2, 2.286379e-2, 2.300225e-2]
    gfa = [0.000000, 5.420000e-8, 0.000000]
    gfb = [9.990000e-9, 0.000000, 9.770000e-9]
    if pick_res[0]:
        index = pick_res[1]
        res_E =  [res_E[index]]
        J     =  [J[index]]     
        gn    =  [gn[index]]    
        gg    =  [gg[index]]    
        gfa   =  [gfa[index]]   
        gfb   =  [gfb[index]]
    ap = 0.948 #[barns]
    A = 238
    
    e_bins = np.logspace(0,2,nbins)
    e_freq = np.zeros(nbins-1)
    scat_freq = np.zeros(nbins-1)
    
    captures = 0
    for i in range(neutrons):
        energy = 1000 
        collisions = 0
        while energy > 1:
            e_freq_step = np.histogram(energy,e_bins)[0]
            e_freq += e_freq_step
            mod_xs_micro = 20
            mod_xs_macro = mod_xs_micro*mf_ratio
            elastic_xs_micro = xs_from_res(ap,A,res_E,J,gn,gg,gfa,gfb,temp,energy,reaction='elastic')
            elastic_xs_macro = elastic_xs_micro
            capture_xs_micro = xs_from_res(ap,A,res_E,J,gn,gg,gfa,gfb,temp,energy,reaction='capture')
            capture_xs_macro = capture_xs_micro
            fuel_xs_macro = elastic_xs_macro+capture_xs_macro
            total_xs_macro = mod_xs_macro + fuel_xs_macro
            mf_random = np.random.uniform(0,1)
            if mf_random < fuel_xs_macro/total_xs_macro:                         #Fuel reaction
                collision_random = np.random.uniform(0,1)
                if collision_random < capture_xs_macro/fuel_xs_macro:     #Capture in fuel
                    captures += 1
                    break
                scat_freq += e_freq_step
                alpha = ((A-1)/(A+1))**2
                energy = (energy-alpha*energy)*np.random.uniform(0,1)+alpha*energy    #Scatter in fuel
            elif mf_random > fuel_xs_macro/total_xs_macro:
                energy = np.random.uniform(0,1)*energy #Scatter in moderator
    e_freq = np.concatenate([[0],e_freq])/neutrons
    scat_freq = np.concatenate([[0],scat_freq])/neutrons
    abs_ratio = captures/neutrons
    if pick_res[0]:
        return e_bins,e_freq,abs_ratio,scat_freq
    else:
        return e_bins,e_freq,abs_ratio

# Plots!
plot_a = False
plot_b = True
plot_c = False
if plot_a:
    ax = plt.subplot(1,1,1) 
    ax.set_xscale("log", nonposx='clip')
    legend_string = []
    neutrons = 200000; logemin = 0; logemax = 2; nbins = 75; temp = 0;
    mf_ratio = 10; 
    e_bins,e_freq,abs_ratio = zerod_mc(mf_ratio, temp, neutrons, logemin, logemax, nbins)
    plt.step(e_bins,e_freq)
    legend_string.append('Mod/Fuel Ratio: '+str(mf_ratio)+', Normalized Absorption='+str(abs_ratio))
    mf_ratio = 1000;
    e_bins,e_freq,abs_ratio = zerod_mc(mf_ratio, temp, neutrons, logemin, logemax, nbins)
    plt.step(e_bins,e_freq)
    legend_string.append('Mod/Fuel Ratio: '+str(mf_ratio)+', Normalized Absorption='+str(abs_ratio))
    mf_ratio = 1000000;
    e_bins,e_freq,abs_ratio = zerod_mc(mf_ratio, temp, neutrons, logemin, logemax, nbins)
    plt.step(e_bins,e_freq)
    legend_string.append('Mod/Fuel Ratio: '+str(mf_ratio)+', Normalized Absorption='+str(abs_ratio))
    plt.xlabel('Energy (eV)')
    plt.ylabel('Normalized Flux (#/src neutron)')
    plt.legend(legend_string)
    plt.title('T='+str(temp)+' K')
    plt.show()

if plot_b:
    ax = plt.subplot(1,1,1)
    ax.set_xscale("log", nonposx='clip')
    legend_string = []
    neutrons = 200000; logemin = 0; logemax = 2; nbins = 75; temp = 1000;
    mf_ratio = 10; 
    e_bins,e_freq,abs_ratio = zerod_mc(mf_ratio, temp, neutrons, logemin, logemax, nbins)
    plt.step(e_bins,e_freq)
    legend_string.append('Mod/Fuel Ratio: '+str(mf_ratio)+', Normalized Absorption='+str(abs_ratio))
    mf_ratio = 1000;
    e_bins,e_freq,abs_ratio = zerod_mc(mf_ratio, temp, neutrons, logemin, logemax, nbins)
    plt.step(e_bins,e_freq)
    legend_string.append('Mod/Fuel Ratio: '+str(mf_ratio)+', Normalized Absorption='+str(abs_ratio))
    mf_ratio = 1000000;
    e_bins,e_freq,abs_ratio = zerod_mc(mf_ratio, temp, neutrons, logemin, logemax, nbins)
    plt.step(e_bins,e_freq)
    legend_string.append('Mod/Fuel Ratio: '+str(mf_ratio)+', Normalized Absorption='+str(abs_ratio))
    plt.xlabel('Energy (eV)')
    plt.ylabel('Normalized Flux (#/src neutron)')
    plt.legend(legend_string)
    plt.title('T='+str(temp)+' K')
    plt.show()

if plot_c:
    ax = plt.subplot(1,1,1)
    ax.set_xscale("log", nonposx='clip')
    legend_string = []
    neutrons = 100000; logemin = 0; logemax = 2; nbins = 75; temp = 1000;
    mf_ratio = 1000;
    res_Es = [6.67, 20.87, 36.68]
    for i in range(3):
        e_bins,e_freq,abs_ratio,scat_freq = zerod_mc(mf_ratio, temp, neutrons, logemin, logemax, nbins, pick_res=[True,i])
        plt.step(e_bins,scat_freq)
        legend_string.append('Resonance at '+str(res_Es[i])+' eV')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Scatters (#/src neutron)')
    plt.legend(legend_string)
    plt.title('T='+str(temp)+' K  M/F Ratio='+str(mf_ratio))
    plt.show()
# e_bin_width = np.zeros([len(e_bins)])
# for i in range(len(e_bins)-2):
#     e_bin_width[i+1] = e_bins[i+2]-e_bins[i+1]
# 
# ax.set_xscale("log", nonposx='clip')
# plt.step(e_bins,e_freq/e_bin_width)
# plt.ylabel('Normalized Flux/eV')
# plt.show()







