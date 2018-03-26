# Script for HW4 Problem 1
# Runs a single pincell to find k_inf with varying pitch
import openmc
import matplotlib.pyplot as plt
import numpy as np
import subprocess

### How to run bash lines from python script
# subprocess.Popen("command here")


#############################
###       MATERIALS       ###
#############################
uo2 = openmc.Material(1, "uo2")
uo2.add_element('U', 1.0, enrichment=3.0)
uo2.add_element('O', 2.0)
uo2.set_density('g/cm3', 10.0)

zirconium = openmc.Material(2, "zirconium")
zirconium.add_element('Zr', 1.0)
zirconium.set_density('g/cm3', 6.6)

water = openmc.Material(3, "h2o")
water.add_nuclide('H1', 2.0)
water.add_element('O', 1.0)
water.set_density('g/cm3', 1.0)
water.add_s_alpha_beta('c_H_in_H2O')

mats = openmc.Materials([uo2, zirconium, water])
mats.export_to_xml()


#############################
###       GEOMETRY        ###
#############################
universe = openmc.Universe()

fuel_or = openmc.ZCylinder(R=0.39)
clad_ir = openmc.ZCylinder(R=0.40)
clad_or = openmc.ZCylinder(R=0.46)

fuel_region = -fuel_or
gap_region = +fuel_or & -clad_ir
clad_region = +clad_ir & -clad_or

fuel = openmc.Cell(1, 'fuel')
fuel.fill = uo2
fuel.region = fuel_region

gap = openmc.Cell(2, 'air gap')
gap.region = gap_region

clad = openmc.Cell(3, 'clad')
clad.fill = zirconium
clad.region = clad_region

pitch = 1.26
box = openmc.get_rectangular_prism(width=pitch, height=pitch,
                                   boundary_type='reflective')
water_region = box & +clad_or
moderator = openmc.Cell(4, 'moderator')
moderator.fill = water
moderator.region = water_region

root = openmc.Universe(cells=(fuel, gap, clad, moderator))
geom = openmc.Geometry(root)
geom.export_to_xml()

#####################################
###        SOURCE/BATCHES         ###
#####################################
point = openmc.stats.Point((0, 0, 0))
src = openmc.Source(space=point)

settings = openmc.Settings()
settings.source = src
settings.batches = 100
settings.inactive = 10
settings.particles = 1000

settings.export_to_xml()


#############################
###       TALLIES         ###
#############################
cell_filter = openmc.CellFilter(fuel)

t = openmc.Tally(1)
t.filters = [cell_filter]
t.nuclides = ['U235']
t.scores = ['total', 'fission', 'absorption', '(n,gamma)']

tallies = openmc.Tallies([t])
tallies.export_to_xml()

openmc.run()


#############################
###       PLOTTING        ###
#############################
p = openmc.Plot()
p.filename = 'pinplot'
p.width = (pitch, pitch)
p.pixels = (200, 200)
p.color_by = 'material'
p.colors = {uo2: 'yellow', water: 'cyan'}

plots = openmc.Plots([p])
plots.export_to_xml()

openmc.plot_geometry()
subprocess.call(["convert","pinplot.ppm","pinplot.png"])
