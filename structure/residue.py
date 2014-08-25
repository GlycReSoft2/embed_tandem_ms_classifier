from composition import Composition

residue_symbol = {
    'A': 'Ala',
    'R': 'Arg',
    'N': 'Asn',
    'D': 'Asp',
    'C': 'Cys',
    'E': 'Glu',
    'Q': 'Gln',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'L': 'Leu',
    'K': 'Lys',
    'M': 'Met',
    'F': 'Phe',
    'P': 'Pro',
    'S': 'Ser',
    'T': 'Thr',
    'W': 'Trp',
    'Y': 'Tyr',
    'V': 'Val'}

residue_table = {
    'Ala': 'C3H5NO',
    'Arg': 'C6H12N4O1',
    'Asn': 'C4H6N2O2',
    'Asp': 'C4H5N1O3',
    'Cys': 'C3H5N1O1S1',
    'Glu': 'C5H7NO3',
    'Gln': 'C5H8N2O2',
    'Gly': 'C2H3N1O1',
    'His': 'C6H7N3O1',
    'Ile': 'C6H11N1O1',
    'Leu': 'C6H11N1O1',
    'Lys': 'C6H12N2O1',
    'Met': 'C5H9N1O1S1',
    'Phe': 'C9H9N1O1',
    'Pro': 'C5H7N1O1',
    'Ser': 'C3H5N1O2',
    'Thr': 'C4H7N1O2',
    'Trp': 'C11H10N2O1',
    'Tyr': 'C9H9N1O2',
    'Val': 'C5H9N1O1',
    'HexNAc': 'C8H13NO5',
    'Hex': 'C6H10O5',
    'dHex': 'C6H10O4',
    'NeuAc': 'C11H17NO8',
    'NeuGc': 'C11H17NO9',
    'Water': 'H2O'}


class Residue:

    """Basic mass values for peptide sequences"""

    def __init__(self, aa=''):
        self.byName(aa)
        #global residue_table
        #self.compo = residue_table.get(aa)
        #self.name = aa
        # if self.compo != None:
        #    self.mass = Composition(self.compo).mass
        # else:
        #    self.mass = 0.0

    def byName(self, name):
        self.compo = residue_table.get(name)
        self.name = name
        if self.compo != None:
            self.mass = Composition(self.compo).mass
        else:
            self.mass = 0.0

    def bySymbol(self, symbol):
        name = residue_symbol[symbol]
        self.byName(name)
        self.symbol = symbol

    def __repr__(self):
        return self.name
