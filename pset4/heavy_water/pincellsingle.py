# Script for HW4 Problem 1
# Runs a single pincell to find k_inf with varying pitch
import openmc
import matplotlib.pyplot as plt
import numpy as np
import subprocess

settings = openmc.Settings()
# Set high tolerance to allow use of lower temperature xs
settings.temperature['tolerance']=10000
settings.temperature['method'] = 'nearest'

#############################
###       MATERIALS       ###
#############################
dy2o3_conc = 0 #weight percent
uo2 = openmc.Material(1, "uo2")
uo2.add_element('U', 1.0, enrichment=3.0)
uo2.add_element('U', 1.0) #natural
uo2.add_element('O', 2.0)
uo2.add_element('Dy', dy2o3_conc*(162.5/210.5), 'wo')
uo2.add_element('O', dy2o3_conc*(48/210.5), 'wo')
uo2.remove_nuclide('U234')
uo2.set_density('g/cm3', 10.0)
uo2.temperature = 900 #kelvin

heavy_water = openmc.Material(2, "d2o")
heavy_water.add_nuclide('H2', 2.0)
heavy_water.add_element('O', 1.0)
heavy_water.set_density('g/cm3', 1.11)
heavy_water.add_s_alpha_beta('c_D_in_D2O')
heavy_water.temperature = 600 #kelvin

mats = openmc.Materials([uo2, heavy_water])
mats.export_to_xml()

#############################
###       GEOMETRY        ###
#############################
universe = openmc.Universe()
    
fuel_or = openmc.ZCylinder(R=1.0)
fuel_region = -fuel_or
fuel_cell = openmc.Cell(1, 'fuel')
fuel_cell.fill = uo2
fuel_cell.region = fuel_region

pitch = 3.0
box = openmc.get_rectangular_prism(width=pitch, height=pitch,
                                   boundary_type='reflective')
water_region = box & +fuel_or
moderator = openmc.Cell(2, 'moderator')
moderator.fill = heavy_water
moderator.region = water_region

root = openmc.Universe(cells=(fuel_cell, moderator))
geom = openmc.Geometry(root)
geom.export_to_xml()

#####################################
###        SOURCE/BATCHES         ###
#####################################
point = openmc.stats.Point((0, 0, 0))
src = openmc.Source(space=point)

settings.source = src
settings.batches = 100
settings.inactive = 10
settings.particles = 1000

settings.export_to_xml()


#############################
###       TALLIES         ###
#############################
# Instantiate an empty Tallies object
tallies_file = openmc.Tallies()

# K-Eigenvalue (infinity) tallies
fiss_rate = openmc.Tally(name='fiss. rate')
fiss_rate.scores = ['nu-fission']
tallies_file.append(fiss_rate)

abs_rate = openmc.Tally(name='abs. rate')
abs_rate.scores = ['absorption']
tallies_file.append(abs_rate)

# Resonance Escape Probability tallies
therm_abs_rate = openmc.Tally(name='therm. abs. rate')
therm_abs_rate.scores = ['absorption']
therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.625])]
tallies_file.append(therm_abs_rate)

# Thermal Flux Utilization tallies
fuel_therm_abs_rate = openmc.Tally(name='fuel therm. abs. rate')
fuel_therm_abs_rate.scores = ['absorption']
fuel_therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.625]),
                                       openmc.CellFilter([fuel_cell])]
tallies_file.append(fuel_therm_abs_rate)

# Fast Fission Factor tallies
therm_fiss_rate = openmc.Tally(name='therm. fiss. rate')
therm_fiss_rate.scores = ['nu-fission']
therm_fiss_rate.filters = [openmc.EnergyFilter([0., 0.625])]
tallies_file.append(therm_fiss_rate)

tallies_file.export_to_xml()

#############################
###       PLOTTING        ###
#############################
p = openmc.Plot()
p.filename = 'pinplot'
p.width = (pitch, pitch)
p.pixels = (200, 200)
p.color_by = 'material'
p.colors = {uo2: 'yellow', heavy_water: 'cyan'}

plots = openmc.Plots([p])
plots.export_to_xml()

openmc.plot_geometry(output=False)
pngstring = 'pinplot{}.png'.format(str(pitch))
subprocess.call(['convert','pinplot.ppm',pngstring])


#############################
###       EXECUTION       ###
#############################
openmc.run()
sp = openmc.StatePoint('statepoint.{}.h5'.format(settings.batches))
# Collect all the tallies 
fiss_rate = sp.get_tally(name='fiss. rate')
fiss_rate_df = fiss_rate.get_pandas_dataframe()
abs_rate = sp.get_tally(name='abs. rate')
abs_rate_df = abs_rate.get_pandas_dataframe()
therm_abs_rate = sp.get_tally(name='therm. abs. rate')
therm_abs_rate_df = therm_abs_rate.get_pandas_dataframe()
fuel_therm_abs_rate = sp.get_tally(name='fuel therm. abs. rate')
fuel_therm_abs_rate_df = fuel_therm_abs_rate.get_pandas_dataframe()
therm_fiss_rate = sp.get_tally(name='therm. fiss. rate')
therm_fiss_rate_df = therm_fiss_rate.get_pandas_dataframe()


# Compute k-infinity
kinf = fiss_rate / abs_rate
kinf_df = kinf.get_pandas_dataframe()

# Compute resonance escape probability
res_esc = (therm_abs_rate)/(abs_rate)
res_esc_df = res_esc.get_pandas_dataframe()

# Compute fast fission factor
fast_fiss = fiss_rate/therm_fiss_rate
fast_fiss_df = fast_fiss.get_pandas_dataframe()

# Compute thermal flux utilization
therm_util = fuel_therm_abs_rate / therm_abs_rate
therm_util_df = therm_util.get_pandas_dataframe()

# Compute neutrons produced per absorption
eta = therm_fiss_rate / fuel_therm_abs_rate
eta_df = eta.get_pandas_dataframe()

print(kinf_df)
print(res_esc_df)
print(fast_fiss_df)
print(therm_util_df)
print(eta_df)

# How to extract values from dataframes
# print(kinf_df['mean'])
# print(kinf_df['std. dev.'])


