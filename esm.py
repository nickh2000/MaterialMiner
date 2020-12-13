from pymatgen.core.structure import Composition
from pymatgen import MPRester
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.ext.matproj import TaskType 
import numpy as np

class enforced_semi_metals:

    def __init__(self, entry): 
        if isinstance(entry, dict):
            self.mp_dict = entry
            self.material_ID = self.mp_dict['ID']
            self.formula = self.mp_dict['formula']
            
        if isinstance(entry, str):
            with MPRester(API_KEY) as m:
                self.get_sg = m.get_space_group_number(entry).as_dict()
                self.mp_dict = self.get_sg
                self.formula = Composition(entry).as_dict()
                self.material_id = entry
                self.formula = entry

    def is_torus(self):
        torus_sg = [1,2,3,5,6,8,10,12,16,21,22,23,25,35,38,42,44,47,65,69,71,75,79,81,82,83,87,89,97,99,107,111,115,119,121,123,139,143,146,147,148,149,150,155,156,157,160,162,164,166,168,174,175,177,183,187,189,191,195,196,197,200,202,204,207,209,211,215,216,217,221,225,229]
        with MPRester() as m:
            
            sg = m.get_space_group_number(self.material_id)
          
            electrons = Composition(self.composition).total_electrons
            if sg in (torus_sg):
                check = electrons % 2
                if (check % 2 == 0): 
                    return m.query(criteria={'task_id':self.material_id}, properties=['pretty_formula'])[0]['pretty_formula']
                    
            else:
                pass 

    def dicosm(self):
        dicosm_sg = [4, 11, 14, 17, 18,20,26,32,36,51,53,55,58,59,63,64,90,94,113,114,127,128,129,135,136,137,173,176,182,185,186,193,194]
        with MPRester() as m:
            sg = m.get_space_group_number(self.material_id)
            electrons = Composition(self.composition).total_electrons
            if sg in (dicosm_sg):
                check = electrons % 8
                if (check % 2 == 0): 
                    return m.query(criteria={'task_id':self.material_id}, properties=['pretty_formula']) 

            else:
                pass 

    def didicosm_1(self): 
        didicosm_sg_1 = [19,61,62,92,198,205,212,213]
        with MPRester() as m:
            sg = m.get_space_group_number(self.material_id)
            electrons = Composition(self.composition).total_electrons
            if sg in (didicosm_sg_1):
                check = electrons % 8
                if (check % 2 == 0): 
                    return m.query(criteria={'task_id':self.material_id}, properties=['pretty_formula'])
            else:
                pass

    def didicosm(self):
        didicosm_sg_2 = [24,73,74,98,122,141,142,199,206,210,214,220,227,228,230]
        with MPRester() as m:
            sg = m.get_space_group_number(self.material_id)
            electrons = Composition(self.composition).total_electrons
            if sg in (didicosm_sg_2):
                 check = electrons % 4
                 if (check % 2 == 0): 
                     return m.query(criteria={'task_id':self.material_id}, properties=['pretty_formula'])
            else:
                pass
    def Amphicosm_1(self):
        amp_sg_1 = [7,13,14,26,27,28,30,31,32,34,39,41,48,49,50,51,53,55,58,59,64,67,68,85,86,100,101,102,103,104,106,116,117,118,124,125,126,127,128,129,132,133,134,135,136,137,201,222,224]
        with MPRester() as m:
            sg = m.get_space_group_number(self.material_id)
            electrons = Composition(self.composition).total_electrons
            if sg in (amp_sg_1):
                 check = electrons % 4
                 if (check % 2 == 0): 
                     return m.query(criteria={'task_id':self.material_id}, properties=['pretty_formula'])

            else:
                pass
    def Amphicosm_2(self):
        amp_sg_2 = [9,15,36,37,40,41,43,45,46,63,64,66,68,70,72,73,74,88,103,104,105,106,108,109,110,112,114,120,122,124,126,128,131,133,135,137,140,141,142,158,159,161,163,165,167,184,185,186,188,190,192,193,194,203,206,218,219,220,222,223,226,227,228,230]
        with MPRester() as m:
            sg = m.get_space_group_number(self.material_id)
            electrons = Composition(self.composition).total_electrons
            if sg in (amp_sg_2):
                 check = electrons % 4
                 if (check % 2 == 0): 
                     return m.query(criteria={'task_id':self.material_id}, properties=['pretty_formula'])[0]['pretty_formula']

            else:
                pass
    def Tetracosm_1(self):
        tet_sg_1 = [76,78,91,92,95,96,212,213]
        with MPRester() as m:
            sg = m.get_space_group_number(self.material_id)
            electrons = Composition(self.composition).total_electrons
            if sg in (tet_sg_1):
                 check = electrons % 8
                 if (check % 2 == 0): 
                     return m.query(criteria={'task_id':self.material_id}, properties=['pretty_formula'])

            else:
                pass 
    def Tetracosm_2(self):
        tet_sg_2 = [77,84,86,93,94,101,102,105,106,131,132,133,134,135,136,137.208,223,224]
        with MPRester() as m:
            sg = m.get_space_group_number(self.material_id)
            electrons = Composition(self.composition).total_electrons
            if sg in (tet_sg_2):
                check = electrons % 4
                if (check % 2 == 0): 
                    return m.query(criteria={'task_id':self.material_id}, properties=['pretty_formula'])
            else:
                pass 

    def Tetracosm_3(self):
        tet_sg_3 = [80,88,98,109,110,141,142,210,214,227,228,230]
        with MPRester() as m:
            sg = m.get_space_group_number(self.material_id)
            electrons = Composition(self.composition).total_electrons
            if sg in (tet_sg_3):
                 check = electrons % 4
                 if (check % 2 == 0): 
                     return m.query(criteria={'task_id':self.material_id}, properties=['pretty_formula'])
            else:
                pass
    def Tricosm(self):
        tri_sg = [144,145,151,152,153,154,171,172,180,181]
        with MPRester() as m:
            sg = m.get_space_group_number(self.material_id)
            electrons = Composition(self.composition).total_electrons
            if sg in (tri_sg):
                 check = electrons % 6
                 if (check % 2 == 0): 
                     return m.query(criteria={'task_id':self.material_id}, properties=['pretty_formula'])
            else:
                pass
    def Hexacosm_1(self):
        hex_sg_1 = [169,170,178,179]
        with MPRester() as m:
            sg = m.get_space_group_number(self.material_id)
            electrons = Composition(self.composition).total_electrons
            if sg in (hex_sg_1):
                 check = electrons % 12
                 if (check % 2 == 0): 
                    return m.query(criteria={'task_id':self.material_id}, properties=['pretty_formula'])
            else:
                pass

    def Hexacosm_2(self):
        hex_sg_2 = [171,172,180,181]
        with MPRester() as m:
            sg = m.get_space_group_number(self.material_id)
            electrons = Composition(self.composition).total_electrons
            if sg in (hex_sg_2):
                 check = electrons % 6
                 if (check % 2 == 0):
                     return m.query(criteria={'task_id':self.material_id}, properties=['pretty_formula'])
            else:
                pass

    def Hexacosm_3(self):
        hex_sg_3 = [173,176,182,185,186,193,194]
        with MPRester() as m:
            sg = m.get_space_group_number(self.material_id)
            electrons = Composition(self.composition).total_electrons
            if sg in (hex_sg_3):
                 check = electrons % 4
                 if (check % 2 == 0): 
                     return m.query(criteria={'task_id':self.material_id}, properties=['pretty_formula'])
            else:
                pass

    def Amphidicosm_1(self):
        amph_sg_1 = [29,54,57,60,61,205]
        with MPRester() as m:
            sg = m.get_space_group_number(self.material_id)
            electrons = Composition(self.composition).total_electrons
            if sg in (amph_sg_1):
                 check = electrons % 8
                 if (check % 2 == 0): 
                     return m.query(criteria={'task_id':self.material_id}, properties=['pretty_formula'])
            else:
                pass
    def Amphicosm_2(self):
        amph_sg_2 = [33,52,56,60,62,130,138]
        with MPRester() as m:
            sg = m.get_space_group_number(self.material_id)
            electrons = Composition(self.composition).total_electrons
            if sg in (amph_sg_2):
                 check = electrons % 8
                 if (check % 2 == 0): 
                     return m.query(criteria={'task_id':self.material_id}, properties=['pretty_formula']) 
            else:
                pass









                    




    






