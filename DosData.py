from pymatgen import MPRester, Spin, Element
from pymatgen import *
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

TiCoO3 = 'mp-19424'
Fe = 'mp-13'

API_KEY = 'UgRqoHkuZyJEVX2d'

max_width = 2
max_bound = 5
precision = .015



'''Class containing information for a material's density of states'''



class DosData:
    
    def __init__(self, arg):

        if isinstance(arg, pd.Series):
            self.efermi = arg['efermi']
            self.energies = np.array(arg['Energies'][1:-1].split(', ')).astype(np.float)
            self.energies = np.subtract(self.energies, self.efermi)
            self.densities = np.array(arg['Densities'][1:-1].split(', ')).astype(np.float)
            self.material_id = arg['ID']

        if isinstance(arg, dict):
            self.dos_dict = arg
            self.material_ID = self.dos_dict['ID']
            self.densities = self.dos_dict['densities']
        if isinstance(arg, str):
            with MPRester(API_KEY) as m:
                self.dos_dict = m.get_dos_by_material_id(arg).as_dict()
                self.densities = self.dos_dict['densities']['1']
                self.efermi = self.dos_dict['efermi']
                self.energies = np.subtract(self.dos_dict['energies'], self.efermi)
                self.material_id = arg
        
        self.dos_array = {energy : density for energy, density in zip(self.energies, self.densities)}


    def get_formula(self):
        with MPRester(API_KEY) as m:
            #print the formula and bandplot for given candidate materials
            return m.query(criteria={'task_id':self.material_id}, properties=['pretty_formula'])[0]['pretty_formula']

    def display_dos_data(self):
        print("The Fermi energy is {}".format(self.dos_obj.efermi))

        par = ["energy", "total_dos", "up_dos", "down_dos"] # arranges a table format for energy, total-dos, spin-up/down dos as the main parameters. 
        
        up_dos = self.dos_obj.get_densities(Spin.up) #define spin-up dos
        down_dos = self.dos_obj.get_densities(Spin.down) #define spin-down dos

        element_dos = self.dos_obj.get_element_dos()   #   defines the total-dos for the specific elements involved
        element_symbols = [e.symbol for e in element_dos]
        par.extend(["{}_dos".format(e) for e in element_symbols]) 

        orbital_dos = self.dos_obj.get_spd_dos() #  defines orbital dos parameter
        orbital_symbols = [o for o in orbital_dos]
        par.extend(["{}_dos".format(o) for o in orbital_symbols]) #adds the orbital-dos to the array for all correspoding energies 

        print ('#'.join(par))  # prints the header for dos data

        for i, energy in enumerate(self.energies): 
            data = [energy, self.dos_array[energy], up_dos[i], down_dos[i]]  #specifies the complete dos data in the new array 
            for e in element_symbols:
                data.append(element_dos[Element(e)].get_densities()[i])   #adds element-dos parameters to the array
            for o in orbital_symbols:
                data.append(orbital_dos[o].get_densities()[i])  #adds orbital-dos parameters to the array

            print ('#'.join([str(d) for d in data])) # prints the complete dos data-set for all corresponing energies 

    #outputs a density of states for energies within the defined bounds
    def center_dos(self, upper_band, lower_band):

        start_energy_index = find_nearest_index(self.energies, lower_band, False)
        end_energy_index = find_nearest_index(self.energies, upper_band, True)

        return {energy: self.dos_array[energy] for energy in self.energies[start_energy_index:end_energy_index]}


    #determines whether the given DOS is sufficiently isolated
    def has_isolated_band(self, upper_band, lower_band, width, offset):

        dos = self.center_dos(upper_band + offset, lower_band + offset)

        for energy, density in dos.items():
            #if the material has densities above or below the spagethii-band widthl, it is not a candidate
            if density != 0 and abs(energy - offset) > width:
                return False

        dosAtWidth = self.center_dos(.1 + offset, -.1 + offset)

        #if theres an occupied energy state within the spaghetti width, return True
        for energy, density in dosAtWidth.items():
            if density != 0:
                return True

        return False

    def get_kpoints(self):
        with MPRester() as m:
            structure = m.get_bandstructure_by_material_id(self.material_id)

        #get bandstructure array for material
        return structure.bands[Spin.up]


    def get_bounds(self):

        fermi_index = find_nearest_index(self.energies, 0, 1)

        max_index = find_nearest_index(self.energies, max_bound, 0)

        min_index = find_nearest_index(self.energies, -max_bound, 2)
        
        upper_width = max_bound
        lower_width = -max_bound
        upper_width_index = max_index
        lower_width_index = min_index
        lower_gap = 0
        upper_gap = 0

        for i in range(fermi_index, max_index):
            if self.densities[i] == 0:
                upper_width = ((int)(self.energies[i - 1] / precision)) * precision
                upper_width_index = i - 1
                break

        for i in np.arange(fermi_index - 1, min_index, -1):
            if self.densities[i] == 0:
                lower_width = ((int)(self.energies[i] / precision)) * precision
                lower_width_index = i
                break

        if self.densities[fermi_index] == 0:
            upper_width = 0
            upper_gap = max_bound
        if self.densities[fermi_index - 1] == 0:
            lower_width = 0
            lower_gap = -max_bound

        #handle the gaps
        for energy in self.energies[upper_width_index + 1 : max_index]:
            upper_gap = round(((int)(energy / precision)) * precision, 4)
            if self.dos_array[energy] != 0:
                break

        for energy in reversed(self.energies[min_index: lower_width_index]):
            lower_gap = round(((int)(energy / precision) - 1 ) * precision, 4)
            if self.dos_array[energy] != 0:
                break

        return upper_width, lower_width, upper_gap, lower_gap


    def plot_dos(self, upper, lower):
        centered_array = self.center_dos(upper, lower)
        plt.plot(list(centered_array.values()), list(centered_array.keys()))
        
        plt.show()

    def get_parameters(self):
        
        bounds = np.array(self.get_bounds())
        
        max_width = bounds[np.argmax(np.abs(bounds[:2]))]

        min_gap = bounds[np.argmin(np.abs(bounds[2:])) + 2]

        return np.abs([max_width, min_gap])
    def get_magmom(self):
        with MPRester(API_KEY) as m:
            structure = mpr.get_structure_by_material_id(self.material_id)
            magmoms = structure.site_properties['magmom'] #compiles all the magnetic momnets into an array 
        return magmoms #return the magmom array 
    
    def verify_magnetic(self): 
        
        return any(map(abs, self.structure.site_properties["magmom"])) #returns true if magnetic moment is present 
    
    def get_nonmagnetic_structure(self, make_primitive: bool = True): 
        structure = self.structure.copy()
        structure.remove_site_property("magmom")
        if make_primitive:
            structure = structure.get_primitive_structure(self.material_id)
            return structure #returns a non-magnetic material



#https://stackoverflow.com/questions/15579649/python-dict-to-numpy-structured-array
def find_nearest_index(array, value, round_up = 0):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    
    #closest value less than desired
    if round_up == 1 and array[idx] < value :
        return idx + 1
    elif round_up == 2 and array[idx] > value:
        return idx - 1
    else:
        return idx
    
