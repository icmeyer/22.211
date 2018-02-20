import numpy as np
import matplotlib.pyplot as plt
from math import log
from res_construction import xs_from_res
from nslowing_resonances import zerod_mc

plot_a = False
plot_b = False
plot_c = True
if plot_a:
    ax = plt.subplot(1,1,1) 
    ax.set_xscale("log", nonposx='clip')
    legend_string = []
    neutrons = 2000; logemin = 0; logemax = 2; nbins = 75; temp = 0;
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
    plt.ylabel('Normalized Flux (#/src neutron/cm^2)')
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
    plt.ylabel('Normalized Flux (#/src neutron/cm^2)')
    plt.legend(legend_string)
    plt.title('T='+str(temp)+' K')
    plt.show()

if plot_c:
    ax = plt.subplot(1,1,1)
    ax.set_xscale("log", nonposx='clip')
    ax.set_yscale("log", nonposx='clip')
    legend_string = []
    neutrons = 100000; logemin = 0; logemax = 2; nbins = 75; temp = 1000;
    mf_ratio = 1000;
    components = ['psi', 'chi', 'pot']
    for i in range(3):
        e_bins,e_freq,abs_ratio,scat_freq = zerod_mc(mf_ratio, temp, neutrons, logemin, logemax, nbins, pick_comp=[True,components[i]])
        plt.step(e_bins,scat_freq)
        legend_string.append(components[i])
    plt.xlabel('Energy (eV)')
    plt.ylabel('Scatters (#/src neutron/cm^2)')
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

