from pymatgen import MPRester
import numpy as np

mpr = MPRester("enter your api key")
structure = mpr.get_structure_by_material_id('mp-768430')
print(structure)
magmoms = structure.site_properties['magmom'] #compiles all the magnetic momnets into an array 
print(magmoms) #prints the magmom array 

