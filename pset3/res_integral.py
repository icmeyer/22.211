# Spin: 0.0 ]
# Scattering length AP: 0.94800 
# 4*pi*AP**2: 11.2934 barns 
# 
# eV          J         GN         GG        GFA        GFB
# ---------- ---------- ---------- ---------- ---------- ----------
# 6.673491+0 5.000000-1 1.475792-3 2.300000-2 0.000000+0 9.990000-9
# 2.087152+1 5.000000-1 1.009376-2 2.286379-2 5.420000-8 0.000000+0
# 3.668212+1 5.000000-1 3.354568-2 2.300225-2 0.000000+0 9.770000-9

#Reconstruct xs from parameters for SLBW
import numpy as np
from math import pi
from scipy.special import wofz
from matplotlib import pyplot as plt
from res_construction import xs_from_res

#xs_d is the dilution xs
def res_integral_inf(energy_min,energy_max,points,temp,lamb,xs_d):
    res_E = [6.673491e+0, 2.087152e+1, 3.668212e+1]
    J = [5.000000e-1, 5.000000e-1, 5.000000e-1]
    gn = [1.475792e-3, 1.009376e-2, 3.354568e-2]
    gg = [2.300000e-2, 2.286379e-2, 2.300225e-2]
    gfa = [0.000000, 5.420000e-8, 0.000000]
    gfb = [9.990000e-9, 0.000000, 9.770000e-9]
    ap = 0.948 #[barns]
    A = 238
    
    energies=np.linspace(energy_min,energy_max,points)
    binwidth = (energy_max-energy_min)/points
    res_int =  0 
    flux_int = 0 
    for i in range(points):
        energy = energies[i]
        capture_xs=xs_from_res(ap,A,res_E,J,gn,gg,
                                         gfa,gfb,temp,energy,
                                         reaction='capture')
        elastic_xs=xs_from_res(ap,A,res_E,J,gn,gg,
                                         gfa,gfb,temp,energy,
                                         reaction='elastic')
        res_int += capture_xs/energy*binwidth
        flux_int += 1/energy*binwidth
    xs_group = res_int/flux_int
        
    return res_int,xs_group

def res_integral_eff(energy_min,energy_max,points,temp,lamb,xs_d):
    res_E = [6.673491e+0, 2.087152e+1, 3.668212e+1]
    J = [5.000000e-1, 5.000000e-1, 5.000000e-1]
    gn = [1.475792e-3, 1.009376e-2, 3.354568e-2]
    gg = [2.300000e-2, 2.286379e-2, 2.300225e-2]
    gfa = [0.000000, 5.420000e-8, 0.000000]
    gfb = [9.990000e-9, 0.000000, 9.770000e-9]
    ap = 0.948 #[barns]
    A = 238
    
    energies=np.linspace(energy_min,energy_max,points)
    binwidth = (energy_max-energy_min)/points
    res_int =  0 
    flux_int = 0 
    for i in range(points):
        energy = energies[i]
        capture_xs=xs_from_res(ap,A,res_E,J,gn,gg,
                                         gfa,gfb,temp,energy,
                                         reaction='capture')
        elastic_xs=xs_from_res(ap,A,res_E,J,gn,gg,
                                         gfa,gfb,temp,energy,
                                         reaction='elastic')
        res_int += capture_xs*(lamb*elastic_xs+xs_d)/(capture_xs+lamb*elastic_xs+xs_d)/energy*binwidth
        flux_int += (lamb*elastic_xs+xs_d)/(capture_xs+lamb*elastic_xs+xs_d)/energy*binwidth
    xs_group = res_int/flux_int
        
    return res_int,xs_group

print("For Problem 3 Infinite Dilute")
print("300 K")
print(res_integral_inf(6,10,10000,300,0,0))
print(res_integral_inf(10,25,10000,300,0,0))
print(res_integral_inf(25,50,10000,300,0,0))
print("1000 K")
print(res_integral_inf(6,10,10000, 1000,0,0))
print(res_integral_inf(10,25,10000,1000,0,0))
print(res_integral_inf(25,50,10000,1000,0,0))
print("For Problem 4: Narrow and Wide models")
print("-------------2000 barns bg---------------")
print("NR lambda = 1")
print(res_integral_eff(6,10,10000,300,1,2000))
print(res_integral_eff(10,25,10000,300,1,2000))
print(res_integral_eff(25,50,10000,300,1,20000))
print("WR lambda = 0")
print(res_integral_eff(6,10,10000,300,0,2000))
print(res_integral_eff(10,25,10000,300,0,2000))
print(res_integral_eff(25,50,10000,300,0,20000))
print("lambda = 0.5")
print(res_integral_eff(6,10,10000,300,0.5,2000))
print(res_integral_eff(10,25,10000,300,0.5,2000))
print(res_integral_eff(25,50,10000,300,0.5,20000))
print("-------------200 barns bg---------------")
print("NR lambda = 1")
print(res_integral_eff(6,10,10000,300,1,200))
print(res_integral_eff(10,25,10000,300,1,200))
print(res_integral_eff(25,50,10000,300,1,200))
print("WR lambda = 0")
print(res_integral_eff(6,10,10000,300,0,200))
print(res_integral_eff(10,25,10000,300,0,200))
print(res_integral_eff(25,50,10000,300,0,200))
print("lambda = 0.5")
print(res_integral_eff(6,10,10000,300,0.5,200))
print(res_integral_eff(10,25,10000,300,0.5,200))
print(res_integral_eff(25,50,10000,300,0.5,200))
print("-------------20 barns bg---------------")
print("NR lambda = 1")
print(res_integral_eff(6,10,10000,300,1,20))
print(res_integral_eff(10,25,10000,300,1,20))
print(res_integral_eff(25,50,10000,300,1,20))
print("WR lambda = 0")
print(res_integral_eff(6,10,10000,300,0,20))
print(res_integral_eff(10,25,10000,300,0,20))
print(res_integral_eff(25,50,10000,300,0,20))
print("lambda = 0.5")
print(res_integral_eff(6,10,10000,300,0.5,20))
print(res_integral_eff(10,25,10000,300,0.5,20))
print(res_integral_eff(25,50,10000,300,0.5,20))
