from pymatgen import MPRester, Spin, Element
from pymatgen import *
import numpy as np
from matplotlib import pyplot as plt

TiCoO3 = 'mp-19424'
Fe = 'mp-13'

API_KEY = 'UgRqoHkuZyJEVX2d'

max_width = 2
max_bound = 5
precision = .015

class DosData:
    def __init__(self, material_id):
        with MPRester(API_KEY) as m:
            self.dos_obj = m.get_dos_by_material_id(material_id)

        self.material_id = material_id

        self.efermi = self.dos_obj.efermi

        self.energies = np.subtract(self.dos_obj.energies, self.efermi)
        
        self.densities = self.dos_obj.get_densities()
        
        self.dos_array = {energy : density for energy, density in zip(self.energies, self.densities)}

    def get_formula(self):
        with MPRester(API_KEY) as m:
            #print the formula and bandplot for given candidate materials
            return m.query(criteria={'task_id':self.material_id}, properties=['pretty_formula'])[0]['pretty_formula']


    def display_dos_data(self):
        print("The Fermi energy is {}".format(self.dos_obj.efermi))

        par = ["energy", "total_dos", "up_dos", "down_dos"]
        
        up_dos = self.dos_obj.get_densities(Spin.up)
        down_dos = self.dos_obj.get_densities(Spin.down)

        element_dos = self.dos_obj.get_element_dos()
        element_symbols = [e.symbol for e in element_dos]
        par.extend(["{}_dos".format(e) for e in element_symbols])

        orbital_dos = self.dos_obj.get_spd_dos()
        orbital_symbols = [o for o in orbital_dos]
        par.extend(["{}_dos".format(o) for o in orbital_symbols])

        print ('#'.join(par)) 

        for i, energy in enumerate(self.energies):
            data = [energy, self.dos_array[energy], up_dos[i], down_dos[i]]
            for e in element_symbols:
                data.append(element_dos[Element(e)].get_densities()[i])  
            for o in orbital_symbols:
                data.append(orbital_dos[o].get_densities()[i]) 

            print ('#'.join([str(d) for d in data]))

    #outputs a density of states for energies within the defined bounds
    def center_dos(self, upper_band, lower_band):

        start_energy_index = find_nearest_energy(self.energies, lower_band, False)
        end_energy_index = find_nearest_energy(self.energies, upper_band, True)

        return {energy: self.dos_array[energy] for energy in self.energies[start_energy_index:end_energy_index]}


    #determines whether the given DOS is sufficiently isolated
    def is_valid(self, upper_band, lower_band, width):

        dos = self.center_dos(upper_band, lower_band)

        for energy, density in dos.items():
            #if the material has densities far enough away from the fermi level, it is not a candidate
            if density != 0 and abs(energy) > width:
                return False
        return True

    def get_kpoints(self):
        with MPRester() as m:
            structure = m.get_bandstructure_by_material_id(self.material_id)

        #get bandstructure array for material
        return structure.bands[Spin.up]

    def get_width_bounds(self):

        for i in np.arange(0, max_width, precision):
            width_upper = self.energies[find_nearest_energy(self.energies, i)]
            if self.dos_array[width_upper] == 0:
                break;

        for i in np.arange(0, -max_width, -precision):
            width_lower = self.energies[find_nearest_energy(self.energies, i)]
            if [width_lower] == 0:
                break;
        return width_upper, width_lower       

    def get_all_bounds(self):

        upper_width, lower_width = self.get_width_bounds()

        for i in np.arange(upper_width + precision, max_bound, precision):
            upper_gap = self.energies[find_nearest_energy(self.energies, i)]
            if self.dos_array[upper_gap] != 0:
                break
        for i in np.arange(lower_width - precision, -max_bound, -precision):

            lower_gap = self.energies[find_nearest_energy(self.energies, i)]
            if self.dos_array[lower_gap] != 0:
                break
        return upper_width, lower_width, upper_gap, lower_gap 


#https://stackoverflow.com/questions/15579649/python-dict-to-numpy-structured-array
def find_nearest_energy(array, value, round_up = True):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    
    #closest value less than desired
    if round_up and array[idx] < value :
        return idx + 1
    elif not round_up and array[idx] > value:
        return idx - 1
    else:
        return idx


if __name__ == "__main__":

    material_id = "mp-720338" # enter mp key

    dos = DosData(material_id)

    dos.display_dos_data()