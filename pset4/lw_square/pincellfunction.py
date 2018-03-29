# Script for HW4 Problem 1
import openmc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess

# dy2o3_conc in atom percent
def pincellfunction(pitch,enrichment): 

    settings = openmc.Settings()
    # Set high tolerance to allow use of lower temperature xs
    settings.temperature['tolerance']=10000
    settings.temperature['method'] = 'nearest'
    settings.temperature['multipole'] = True
    settings.cutoff={'energy': 1e-8} #energy cutoff in eV
    
    #############################
    ###       MATERIALS       ###
    #############################
    uo2 = openmc.Material(1, "uo2")
    uo2.add_element('U', 1.0, enrichment=enrichment)
    uo2.add_element('O', 2.0)
    uo2.remove_nuclide('U234')
    uo2.set_density('g/cm3', 10.97)
    uo2.temperature = 900 #kelvin
    
    water = openmc.Material(3, "h2o")
    water.add_nuclide('H1', 2.0)
    water.add_nuclide('O16', 1.0)
    water.set_density('g/cm3', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')
    water.temperature = 600 #kelvin
    
    mats = openmc.Materials([uo2, water])
    mats.export_to_xml()
    
    #############################
    ###       GEOMETRY        ###
    #############################
    universe = openmc.Universe()
        
    fuel_or = openmc.ZCylinder(R=0.5)
    fuel_region = -fuel_or
    fuel_cell = openmc.Cell(1, 'fuel')
    fuel_cell.fill = uo2
    fuel_cell.region = fuel_region
    
    box = openmc.get_rectangular_prism(width=pitch, height=pitch,
                                       boundary_type='reflective')
    water_region = box & +fuel_or
    moderator = openmc.Cell(2, 'moderator')
    moderator.fill = water
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
    p.colors = {uo2: 'yellow', water: 'blue'}
    
    plots = openmc.Plots([p])
    plots.export_to_xml()
    
    openmc.plot_geometry(output = False)
    pngstring = 'pinplot{}.png'.format(str(pitch))
    subprocess.call(['convert','pinplot.ppm',pngstring])
    subprocess.call(['mv',pngstring,'figures/'+pngstring])
    
    
    #############################
    ###       EXECUTION       ###
    #############################
    openmc.run(output=False)
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
    

    columns = ['pitch','enrichment',
             'kinf mean','kinf sd','res_esc mean','res_esc sd',
             'fast_fiss mean','fast_fiss sd','therm_util mean',
             'therm_util sd','eta mean','eta sd']
    data = [[pitch,enrichment,
            kinf_df['mean'][0],kinf_df['std. dev.'][0],
            res_esc_df['mean'][0],res_esc_df['std. dev.'][0],
            fast_fiss_df['mean'][0],fast_fiss_df['std. dev.'][0],
            therm_util_df['mean'][0],therm_util_df['std. dev.'][0],
            eta_df['mean'][0],eta_df['std. dev.'][0]]]
    all_tallies = pd.DataFrame(data, columns = columns)
    
    return all_tallies

    
    # How to extract values from dataframes
    # print(kinf_df['mean'])
    # print(kinf_df['std. dev.'])


