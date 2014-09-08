import re

component = {
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


class Composition:

    """Take composition formula and calculate corresponding mass"""

    def __init__(self, str):
        # Needs update to support general expression
        # This pattern will match an uppercase character followed by 0 or more
        # lowercase characters denoting the element/group, followed by 0 or more
        # digits denoting the quantity of the entity (defaults to 1 if no number is present)
        compo_pattern = re.compile(r'((?:[A-Z][a-z]*)|e|p)(\d*)')

        self.dict = {}
        self.mass = 0.0
        self.compo = str

        for match in compo_pattern.findall(str):
            count = 0
            if match[1] == '':
                count = 1
            else:
                count = int(match[1])
            if match[0] in self.dict:
                self.dict[match[0]] += count
            else:
                self.dict[match[0]] = count
            self.mass += component[match[0]] * count
