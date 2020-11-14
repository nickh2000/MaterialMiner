import numpy as np
from pymatgen import MPRester
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer

if __name__ == "__main__":
    MAPI_KEY = "" 

mpr = MPRester(MAPI_KEY)
s = mpr.get_structure_by_material_id("mp-775157")
spa = SpacegroupAnalyzer(s, symprec=0.1) 
print(spa.get_crystal_system()) 
print(spa.get_point_group_symbol())
print(spa.get_space_group_number())
print(spa.get_space_group_symbol())