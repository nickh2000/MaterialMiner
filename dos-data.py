from pymatgen import MPRester, Spin, Element
from pymatgen import *
import numpy as np
from matplotlib import pyplot as plt


def dos_data(dos):
    print("The Fermi energy is {}".format(dos.efermi))

    par = ["energy", "total_dos", "up_dos", "down_dos"]
    energies = dos.energies
    total_dos = dos.get_densities()  
    up_dos = dos.get_densities(Spin.up)
    down_dos = dos.get_densities(Spin.down)

    element_dos = dos.get_element_dos()
    element_symbols = [e.symbol for e in element_dos]
    par.extend(["{}_dos".format(e) for e in element_symbols])

    orbital_dos = dos.get_spd_dos()
    orbital_symbols = [o for o in orbital_dos]
    par.extend(["{}_dos".format(o) for o in orbital_symbols])

    print ('#'.join(par)) 

    for i, energy in enumerate(energies):
        data = [energy, total_dos[i], up_dos[i], down_dos[i]]
        for e in element_symbols:
            data.append(element_dos[Element(e)].get_densities()[i])  
        for o in orbital_symbols:
            data.append(orbital_dos[o].get_densities()[i]) 

        print ('#'.join([str(d) for d in data]))


if __name__ == "__main__":
  
    API_KEY = "NONE"  # enter API key 
    material_id = "mp-720338" # enter mp key

    mpr = MPRester(API_KEY)
    dos = mpr.get_dos_by_material_id(material_id)

    dos_data(dos)