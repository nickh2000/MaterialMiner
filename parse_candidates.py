from pymatgen import *
import random
from DosData import *
import requests
import pandas as pd
import itertools
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import periodictable 
import json
from bokeh.io import output_file, show
from bokeh.plotting import figure
from bokeh.sampledata.periodic_table import elements
from bokeh.transform import dodge, linear_cmap, log_cmap
import bokeh.palettes
from bokeh.models import ColorBar
import math
from pprint import pprint


API_KEY = '9c8rlQwkCJdE56MD'
mpr = MPRester(API_KEY)

def analyze_candidates(ID):
	ID = ID.rstrip()
	DOS = DosData(ID)			
	centered_array  = DOS.center_dos(10, -10)
	plt.plot(list(centered_array.values()), list(centered_array.keys()))
	plt.ylabel(DOS.get_formula())
	plt.show()

# candidates_select() returns a list of candidates including/excluding specific elements.
# We exclude O and F, 3d block transition metals, f elements 
def candidates_select():
    with open("candidates.txt") as f:
        materials = [material.rstrip() for material in f.readlines()]
        response = mpr.query(criteria={'task_id': {'$in': list(materials)},"elements": {"$in":["B","Si","Ge","Sn","Bi","Mo","Ru","Rh","Ba"], "$nin":["O","F"]}},#, "nelements": 2
        properties=["pretty_formula", "task_id"])         
        pprint(response)

if __name__ == "__main__":
    candidates_select()

# 1) Naive run with common elements found in topological materials: "Sb","Bi","Sn","Te","Pb"

# 2) PHYSICAL REVIEW B 101, 245117 (2020) Elements commonly present in TIs:  
# "B" 17%, "Si" 18%, "Ge" 21%, Sn 16%, Bi 12%, Mo 18%, Ru %20, Rh %20, Ba 15%, 
# ["B","Si","Ge","Sn","Bi","Mo","Ru","Rh","Ba"] 
