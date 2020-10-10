
from pymatgen.electronic_structure.plotter import DosPlotter, BSPlotter
from pymatgen import *
import numpy as np
from matplotlib import pyplot as plt


TiCoO3 = 'mp-19424'
Fe = 'mp-13'

materials_validity = [(TiCoO3, True), (Fe, False), ]


lower_band = -.8
upper_band = 1.2
threshold = .3

test_dos = {-2: 1, -1 : .1, 0:.3, 1: 0}


#https://stackoverflow.com/questions/15579649/python-dict-to-numpy-structured-array
def find_nearest(array, value, round_up = True):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    
    #closest value less than desired
    if round_up and array[idx] < value :
    	return idx + 1
    elif not round_up and array[idx] > value:
    	return idx - 1
    else:
    	return idx

def get_dos_array(material_id):

	with MPRester(api_key="UgRqoHkuZyJEVX2d") as m:

		dos = m.get_dos_by_material_id(TiCoO3)
		energies = dos.as_dict()['energies']
		e_fermi = dos.as_dict()['efermi']
		densities = dos.get_densities()


	return {energy - e_fermi : density for energy, density in zip(energies, densities)}


def center_dos(dos):

	energies = list(dos.keys())

	start_energy_index = find_nearest(energies, lower_band, False)
	end_energy_index = find_nearest(energies, upper_band, True)

	return {energy: dos[energy] for energy in energies[start_energy_index:end_energy_index]}

def is_valid_dos(dos):
	
	for energy, density in dos.items():
		if density != 0 and abs(energy) > threshold:
			return False
	return True

def is_valid_material(material_id):
	return is_valid_dos(get_dos_array(material_id))


if __name__ == "__main__":
	dos = get_dos_array(Fe)
	dos = center_dos(dos)

	print(is_valid_dos(dos))
	plt.plot(list(dos.keys()), list(dos.values()))
	plt.show()
	