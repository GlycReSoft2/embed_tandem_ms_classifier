# Not used?
from composition import Composition

mod_table = {
    'Carbamidomethyl': ['Cys', 57.0214],
    'Deamidated': ['Asn', 0.9840099999999978],
    'HexNAc': ['Asn', 203.07937],
    'pyroGlu': ['Gln', -17.02655]}


class Modification:

    """description of class"""

    def __init__(self, mod_name, mod_pos=-1, mod_num=1, mass=0.0, target=''):
        self.name = mod_name
        self.position = mod_pos
        self.number = mod_num
        key = mod_table.get(mod_name)
        if key == None:
            self.mass = mass
            self.target = target
            mod_table[mod_name] = [target, mass]
        else:
            self.mass = mod_table.get(mod_name)[1]
            self.target = mod_table.get(mod_name)[0]
