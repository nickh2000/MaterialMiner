from pymatgen import MPRester
from pymatgen.electronic_structure.plotter import BSPlotter

if __name__ == "__main__":
    MAPI_KEY = "" #put your api key here as a string
    MP_ID = ""   # enter material id  ex: mp-19017

    mpr = MPRester(MAPI_KEY) 

    my_bs = mpr.get_bandstructure_by_material_id(MP_ID)
    
    plot = "band plot: {}".format(my_bs.get_band_gap())
    print(plot)  
    BSPlotter(my_bs).show()