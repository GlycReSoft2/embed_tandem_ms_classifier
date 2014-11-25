import re
from ..utils.memoize import memoize

# This pattern will match an uppercase character followed by 0 or more
# lowercase characters denoting the element/group, followed by 0 or more
# digits denoting the quantity of the entity (defaults to 1 if no number is present)
composition_formula_pattern = re.compile(r'((?:[A-Z][a-z]*)|e|p)(\d*)')

name_to_mass = {
    'C': 12.00000,
    'H': 1.0078250350,
    'N': 14.0030740000,
    'O': 15.9949146300,
    'S': 31.9720707000,
    'P': 30.9737620000,
    'I': 126.904473,
    'Cl': 34.9688527300,
    'F': 18.9984032200,
    'Br': 78.918361000,
    'Na': 22.989767700,
    'K': 38.9637069000,
    'Ca': 39.9625912000,
    'e': 0.000549,
    'p': 1.00727647
}


@memoize()
def composition_to_mass(formula):
    '''Fast, low allocation mass computation'''
    mass = 0.0
    try:
        tokens = composition_formula_pattern.findall(formula)
    except TypeError, e:
        print("An error ocurred with formula %s" % formula)
        raise e
    try:
        for name, count in tokens:
            if count == '':
                count = 1
            else:
                count = int(count)
            mass += (name_to_mass[name] * count)
    except Exception, e:
        print("An error occurred for name %s, count %s" % (name, count))
        raise
    return mass


class Composition:
    """Take composition formula and calculate corresponding mass"""
    def __init__(self, formula):
        self.dict = {}
        self.mass = 0.0
        self.compo = formula

        for match in composition_formula_pattern.findall(formula):
            count = 0
            if match[1] == '':
                count = 1
            else:
                count = int(match[1])
            if match[0] in self.dict:
                self.dict[match[0]] += count
            else:
                self.dict[match[0]] = count
            self.mass += name_to_mass[match[0]] * count
