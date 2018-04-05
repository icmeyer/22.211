import openmc


pitch =2.0

fuel_or = openmc.ZCylinder(R=0.5)
fuel_region = -fuel_or
fuel_cell = openmc.Cell(1, 'fuel')
fuel_cell.region = fuel_region

hexagon = openmc.get_hexagonal_prism(edge_length = 1/3**(1/2)*pitch,
                                     orientation = 'x',
                                     boundary_type = 'reflective')
bottom = openmc.ZPlane(z0=-pitch/2, boundary_type='reflective')
top = openmc.ZPlane(z0=pitch/2, boundary_type='reflective')

water_region = +bottom & -top & hexagon & +fuel_or
moderator = openmc.Cell(2, 'moderator')
moderator.region = water_region

root = openmc.Universe(cells=(fuel_cell, moderator))
geom = openmc.Geometry(root)
geom.export_to_xml()
print(water_region)

print('(0,0,0)')
print((0,0,0) in fuel_region,(0,0,0) in water_region)
print('(0,0.6,0)')
print((0,0.6,0) in fuel_region,(0,0.6,0) in water_region)
print('(-pitch,pitch,0)')
print((-pitch,pitch,0) in fuel_region,(-pitch,pitch,0) in water_region)

universe = openmc.Universe(cells=(fuel_cell, moderator))
universe = openmc.Universe(cells=[moderator])
universe.plot(width=(2*pitch, 2*pitch))
