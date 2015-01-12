import csv
import os
import re

from collections import OrderedDict

from .modification import AnonymousModificationRule

mammalian_glycomedb_nlinked_path = os.path.join(os.path.dirname(__file__), "data", "Mammalian_GlycomeDB_NLinked.csv")

# http://www.google.com/patents/US20140117225
oxonium_ions = {
    "HexNAc": 204.0864,
    "HexNAc-H2O": 186.0754,
    "HexNAc-2H2O": 168.0650,
    "FragmentOfHexNAc": 138.0542,
    "dHex": 146.05791,
    "Hex": 163.0601,
    "HexHexNAc": 366.1394,
    "NeuAc": 292.1026,
    "Neu5Ac-H2O": 274.0920,
}


class Glycan(object):

    @staticmethod
    def glycan_composition_string(count_list):
        return "[{0}]".format(";".join(map(lambda x: str(int(x)), count_list)))

    def __init__(self, mass, composition, glycan_identities=None, **kwargs):
        if isinstance(composition, basestring):
            composition = map(int, composition[1:-1].split(";"))
        self.composition = composition
        self.mass = mass
        self.glycan_identities = glycan_identities
        self.data = kwargs

        self.composition_dict = OrderedDict()
        if glycan_identities is not None:
            iter_identities = iter(glycan_identities)
            for count in composition:
                self.composition_dict[iter_identities.next()] = count

    @property
    def molecular_weight(self):
        return self.mass

    def as_modification(self):
        return AnonymousModificationRule("Glycan[{0}]".format(';'.join(map(str, self.composition))), self.mass)

    def __repr__(self):
        return "Glycan({composition} {mass})".format(**self.__dict__)

    def __hash__(self):
        return hash(str(self))


class GlycanHypothesis(list):
    def __init__(self, csv_path=None):
        if csv_path is not None:
            self.reader = csv.DictReader(open(csv_path))
            self.glycan_identities = self.get_glycan_identities(self.reader.fieldnames)
            self.extend(self.reader)
        self.extract_glycans()

    def extract_glycans(self):
        self.glycans = list({Glycan(float(r["Molecular Weight"]),
                                    r["Compositions"],
                                    glycan_identities=self.glycan_identities) for r in self})
        return self.glycans

    def get_glycan_identities(self, colnames):
        glycan_identity = []
        extract_state = False
        for col in colnames:
            if col == "Compositions":
                extract_state = True
                continue
            elif col == "Adduct/Replacement":
                extract_state = False
                break
            elif extract_state:
                glycan_identity.append(col.replace("G:", ""))
        return glycan_identity


def from_predictions(frame):
    return frame.apply(lambda x: Glycan(x.glycanMass, x.Glycan), 1).tolist()


def load_from_file(path_to_file=mammalian_glycomedb_nlinked_path):
    gh = GlycanHypothesis(path_to_file)
    return gh.glycans


def pack_glycan_string(glycan_map, glycan_order):
    pack_string = "[{0}]".format(";".join([glycan_map[g] for g in glycan_order]))
    return pack_string


def parse_glycan_pack_string(pack_string, glycan_order):
    glycan_map = {g: int(c) for g, c in zip(glycan_order, re.findall('\\b\\d+\\b', pack_string))}
    return glycan_map
