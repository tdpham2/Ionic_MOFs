import re
from collections import defaultdict
import json

metals = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr', 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra',
          'Al', 'Ga', 'Ge', 'In', 'Sn', 'Sb', 'Tl', 'Pb', 'Bi', 'Po',
          'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
          'Y' , 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
          'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
          'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'U', 'Tm', 'Yb', 'Lu',
          'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'
         ]

mass_key_round = {'X': 0.0, 'NA': 0.0, 'H': 1.01, 'He': 4.0, 'Li': 6.94, 'Be': 9.01, 'B': 10.81, 'C': 12.01, 'N': 14.01, 'O': 16.0, 'F': 19.0, 'Ne': 20.18, 'Na': 22.99, 'Mg': 24.3, 'Al': 26.98, 'Si': 28.09, 'P': 30.97, 'S': 32.06, 'Cl': 35.45, 'Ar': 39.95, 'K': 39.1, 'Ca': 40.08, 'Sc': 44.96, 'Ti': 47.87, 'V': 50.94, 'Cr': 52.0, 'Mn': 54.94, 'Fe': 55.84, 'Co': 58.93, 'Ni': 58.69, 'Cu': 63.55, 'Zn': 65.38, 'Ga': 69.72, 'Ge': 72.64, 'As': 74.92, 'Se': 78.96, 'Br': 79.9, 'Kr': 83.8, 'Rb': 85.47, 'Sr': 87.62, 'Y': 88.91, 'Zr': 91.22, 'Nb': 92.91, 'Mo': 95.96, 'Tc': 98, 'Ru': 101.07, 'Rh': 102.91, 'Pd': 106.42, 'Ag': 107.87, 'Cd': 112.41, 'In': 114.82, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.6, 'I': 126.9, 'Xe': 131.29, 'Cs': 132.91, 'Ba': 137.33, 'La': 138.91, 'Ce': 140.12, 'Pr': 140.91, 'Nd': 144.24, 'Pm': 145, 'Sm': 150.36, 'Eu': 151.96, 'Gd': 157.25, 'Tb': 158.93, 'Dy': 162.5, 'Ho': 164.93, 'Er': 167.26, 'Tm': 168.93, 'Yb': 173.05, 'Lu': 174.97, 'Hf': 178.49, 'Ta': 180.95, 'W': 183.84, 'Re': 186.21, 'Os': 190.23, 'Ir': 192.22, 'Pt': 195.08, 'Au': 196.97, 'Hg': 200.59, 'Tl': 204.38, 'Pb': 207.2, 'Bi': 208.98, 'Po': 209, 'At': 210, 'Rn': 222, 'Fr': 223, 'Ra': 226, 'Ac': 227, 'Th': 232.04, 'Pa': 231.04, 'U': 238.03, 'Np': 237, 'Pu': 244, 'Am': 243, 'Cm': 247, 'Bk': 247, 'Cf': 251, 'Es': 252, 'Fm': 257, 'Md': 258, 'No': 259, 'Lr': 262, 'Rf': 265, 'Db': 268, 'Sg': 271, 'Bh': 272, 'Hs': 270, 'Mt': 276, 'Ds': 281, 'Rg': 280, 'Cn': 285, 'Uut': 284, 'Uuq': 289, 'Uup': 288, 'Uuh': 293, 'Uuo': 294, 'D':2}

class MOFChemicalFormula():
    """ A class for MOF chemical formula obtained from CSD database or ASE. Mostly useful for CSD MOFs

    """
    def __init__(self, formula):
        self.formula = formula

    def format_formula(self, formula, pattern):
        """ Remove a pattern from formula """
        return re.sub(pattern, '', formula)

    def split_formula(self, formula, sep=','):
        """ Split formula by a separator"""
        return formula.split(sep)

    def is_charged(self, formula):
        """ Check if formula is charged by whether a '+' or '-' sign is in the formula. Return True and charge sign if True, False and empty string if False"""
        if '-' in formula:
            return True, '-'
        elif '+' in formula:
            return True, '+'
        else:
            return False, ''

    def count_charge(self, formula):
        """ Return the number of charge of the formula"""
        charged, sign = self.is_charged(formula)

        if charged == True:
            charges = re.findall(r'\d+\{}'.format(sign), formula)
            if len(charges) != 1:
                print("Error! Multiple charges in {}".format(formula))
            else:
                ncharge = float(charges[0][:-1])
                return ncharge, sign
        else:
            return 0, ''

    def element_dict(self, fragment_formula=None):
        """ Return a dictionary of element and its frequency.
        Parameters:
            fragment_formula: str, if None, return the element dict of the whole chemical formula. Else, return the element dict of the input fragment
        """
        element_dict = defaultdict(lambda:0)

        if fragment_formula == None:
            fragment_formula = self.formula

        element_list = re.findall(r"[A-Z]{1}[a-z]{0,1}\d+", fragment_formula)
        for element_item in element_list:
            element_freq = re.match(r"(\D+)(\d+)", element_item)
            element = element_freq.group(1)
            freq = element_freq.group(2)                    
            element_dict[element] += int(freq)

        return element_dict

    def get_mass_from_element_dict(self, element_dict = None):
        """ Return the total mass from element_dict"""

        mass = 0
        if element_dict == None:
            element_dict = self.element_dict()
        for element in element_dict:
            mass += mass_key_round[element] * element_dict[element]
        return round(mass, 2)

    def is_framework(self):
        """ Return the largest fragment that contains metal element"""
        
        return True
    def formula_details(self):
        """ Return details of the formula, including framework formula, solvent and ions"""
        frameworks = {}

        fragments = self.split_formula(self.formula, ',')
        for fragment in fragments:
            #format_fragment = self.format_formula(fragment, '\(|\)')
            format_fragment = fragment
            # Get charge of fragment
            ncharge, sign = self.count_charge(format_fragment)
            
            # Get mass of fragment
            element_dict = self.element_dict(fragment_formula=format_fragment)
            mass = self.get_mass_from_element_dict(element_dict=element_dict)

            frameworks[fragment] = [format_fragment, ncharge, sign, mass]
        return frameworks
    def get_largest_fragment(self, has_metal=True):
        """ Return the largest fragment from a chemical formula
        Parameters: has_metal: boolean, whether the largest fragment must contain a metal element.

        """
        frameworks = self.formula_details()
        largest_mass = 0
        largest_frag = None

        for index, framework_item in enumerate(sorted(list(frameworks))):
            element_dict = self.element_dict(fragment_formula = framework_item)
            for element in element_dict:
                if element in metals:
                    if frameworks[framework_item][-1] > largest_mass:
                        largest_frag = framework_item
                        largest_mass = frameworks[framework_item][-1]
                else:
                    continue

        return largest_frag

    def is_ionic(self):
        """ Check if the framework is ionic or not by determining whether the largest fragment is charged or not"""
        
        frameworks = self.formula_details()
        largest_mass = 0
        largest_frag = None

        framework_dict = {'framework':[],
                          'solvent':[],
                          'cation': [],
                          'anion': []
                         }
        largest_frag = self.get_largest_fragment(has_metal=True)

        if frameworks[largest_frag][2] == '+' or frameworks[largest_frag][2] == '-':
            # Loop again to separate main framework, ion and solvent
            for index, framework_item in enumerate(sorted(list(frameworks))):
                if framework_item == largest_frag:
                    temp_dict = { 'formula': frameworks[framework_item][0],
                                                    'charge': frameworks[framework_item][2],
                                                    'net_charge_per_formula': frameworks[framework_item][1],
                                                    'formula_mass': frameworks[framework_item][3]
                                                    }
 
                    framework_dict['framework'].append(temp_dict)

                else:
                    temp_dict = { 'formula': frameworks[framework_item][0],
                                'charge': frameworks[framework_item][2],
                                'net_charge_per_formula': frameworks[framework_item][1],
                                'formula_mass': frameworks[framework_item][3]
                                }
                    if temp_dict['charge'] == '':
                        framework_dict['solvent'].append(temp_dict)
                    elif temp_dict['charge'] == '-':
                        framework_dict['anion'].append(temp_dict)
                    elif temp_dict['charge'] == '+':
                        framework_dict['cation'].append(temp_dict)
            return True, framework_dict
        else:
            for index, framework_item in enumerate(sorted(list(frameworks))):
                if framework_item == largest_frag:
                    temp_dict = { 'formula': frameworks[framework_item][0],
                                                    'charge': 0,
                                                    'net_charge_per_formula': 0,
                                                    'formula_mass': frameworks[framework_item][3]
                                                    }

                    framework_dict['framework'].append(temp_dict)

                else:
                    temp_dict = { 'formula': frameworks[framework_item][0],
                                'charge': frameworks[framework_item][2],
                                'net_charge_per_formula': frameworks[framework_item][1],
                                'formula_mass': frameworks[framework_item][3]
                                }
                    if temp_dict['charge'] == '':
                        framework_dict['solvent'].append(temp_dict)
                    elif temp_dict['charge'] == '-':
                        framework_dict['anion'].append(temp_dict)
                    elif temp_dict['charge'] == '+':
                        framework_dict['cation'].append(temp_dict)
            return False, framework_dict
    
    def is_multiplier(self, new_element_dict):
        """ Compare two element dictionary to check if new one is a multiplier.
        Parameters:
            new_element_dict: dict, dictionary of the other fragment elements. Same format that is called from element_dict() method.
        """
        element_dict = self.element_dict()
        check = False
        if set(element_dict.keys()) == set(new_element_dict.keys()):
            multipliers = []
            for element in element_dict:
                if new_element_dict[element] % element_dict[element] == 0:
                    multipliers.append(int(new_element_dict[element]/element_dict[element]))
                else:
                    check = False
                    multipliers = []
                    break
            mulipliers_set = set(multipliers)
            if len(mulipliers_set) == 1:
                check = True 
                multiplier = multipliers[0]
            else:
                check =False
                multiplier = 0
        else:
            check = False
            multiplier = 0
        return check, multiplier
