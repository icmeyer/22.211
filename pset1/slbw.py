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


res_E = [6.673491e+0, 2.087152e+1, 3.668212e+1]
J = [5.000000e-1, 5.000000e-1, 5.000000e-1]
gn = [1.475792e-3, 1.009376e-2, 3.354568e-2]
gg = [2.300000e-2, 2.286379e-2, 2.300225e-2]
gfa = [0.000000, 5.420000e-8, 0.000000]
gfb = [9.990000e-9, 0.000000, 9.770000e-9]
ap = 0.948 #[barns]
A = 238

points = 100000
energies = np.linspace(0.01,100,points)
capture_xs = np.zeros([points])
elastic_xs = np.zeros([points])
total_xs = np.zeros([points])
temp = 0
for j in range(points):
    capture_xs[j]=xs_from_res(ap,A,res_E,J,gn,gg,
                                     gfa,gfb,temp,energies[j],
                                     reaction='capture')
    elastic_xs[j]=xs_from_res(ap,A,res_E,J,gn,gg,
                                     gfa,gfb,temp,energies[j],
                                     reaction='elastic')
total_xs =capture_xs+elastic_xs

plt.clf()
plt.cla()
plt.close()

plt.subplot(2,1,1)
plt.loglog(energies,capture_xs,'C0')
plt.loglog(energies,elastic_xs,'C1')
plt.loglog(energies,total_xs,'C2')
leg = plt.legend(['Radiative Capture','Elastic Scattering','Total'])
leg.legendHandles[0].set_color('C0')
leg.legendHandles[1].set_color('C1')
leg.legendHandles[2].set_color('C2')
plt.xlabel('Energy [eV]')
plt.ylabel('Cross Section [barns]')
plt.title('0 K Cross Sections')


capture_xs = np.zeros([points])
elastic_xs = np.zeros([points])
total_xs = np.zeros([points])
temp = 1000
for j in range(points):
    capture_xs[j]=xs_from_res(ap,A,res_E,J,gn,gg,
                                     gfa,gfb,temp,energies[j],
                                     reaction='capture')
    elastic_xs[j]=xs_from_res(ap,A,res_E,J,gn,gg,
                                     gfa,gfb,temp,energies[j],
                                     reaction='elastic')
total_xs = capture_xs+elastic_xs

plt.subplot(2,1,2)
plt.loglog(energies,capture_xs,'C0')
plt.loglog(energies,elastic_xs,'C1')
plt.loglog(energies,total_xs,'C2')
leg = plt.legend(['Radiative Capture','Elastic Scattering','Total'])
leg.legendHandles[0].set_color('C0')
leg.legendHandles[1].set_color('C1')
leg.legendHandles[2].set_color('C2')
plt.xlabel('Energy [eV]')
plt.ylabel('Cross Section [barns]')
plt.title('1000 K Cross Sections')
plt.show()
plt.clf()
plt.cla()
plt.close()

#Plot single cross section to check
points = 100000
energies = np.logspace(-1,3,points)
capture_xs = np.zeros([points])
elastic_xs = np.zeros([points])
total_xs = np.zeros([points])
temp = 0
for j in range(points):
    capture_xs[j]=xs_from_res(ap,A,res_E,J,gn,gg,
                                     gfa,gfb,temp,energies[j],
                                     reaction='capture')
    elastic_xs[j]=xs_from_res(ap,A,res_E,J,gn,gg,
                                     gfa,gfb,temp,energies[j],
                                     reaction='elastic')
total_xs =capture_xs+elastic_xs
plt.loglog(energies,elastic_xs)
plt.show()
print(xs_from_res(ap,A,res_E,J,gn,gg,gfa,gfb,0.0001,6.67,reaction='elastic'))
print(sum(total_xs >=0))
print(sum(capture_xs >=0))
print(sum(elastic_xs >=0))
