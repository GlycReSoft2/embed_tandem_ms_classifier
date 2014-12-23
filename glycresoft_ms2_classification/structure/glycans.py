import csv
import os
import re

from collections import OrderedDict
from operator import itemgetter

from .modification import AnonymousModificationRule

mammalian_glycomedb_nlinked_path = os.path.join(os.path.dirname(__file__), "data", "Mammalian_GlycomeDB_NLinked.csv")


class SimpleGlycan(tuple):
    __slots__ = ()
    _fields = ("molecular_weight", "composition")

    def __new__(cls, molecular_weight, composition):
        return tuple.__new__(cls, (molecular_weight, composition))

    @classmethod
    def _make(cls, iterable, new=tuple.__new__, len=len):
        result = new(cls, iterable)
        if len(result) != 2:
            raise TypeError('Expected 2 arguments, got %d' % len(result))
        return result

    def __repr__(self):
        return "SimpleGlycan(molecular_weight=%r, composition=%r)" % self

    def _asdict(self):
        return OrderedDict(zip(self._fields, self))

    def _replace(self, **kwargs):
        results = self._make(map(kwargs.pop, self._fields, self))
        if kwargs:
            raise ValueError('Got unexpected field names: %r' % kwargs.keys())
        return results

    def __getnewargs__(self):
        'Return self as a plain tuple. Used by copy and pickle.'
        return tuple(self)

    __dict__ = property(_asdict)

    def __getstate__(self):
        'Exclude the OrderedDict from pickling'
        pass

    @property
    def mass(self):
        return self.molecular_weight

    def as_modification(self):
        return AnonymousModificationRule("Glycan" + self.composition, self.molecular_weight)

    molecular_weight = property(itemgetter(0), doc='Alias for field number 0')

    composition = property(itemgetter(1), doc='Alias for field number 1')


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


def glycan_from_predictions(classify_matches_data):
    return classify_matches_data.apply(lambda x: Glycan(x.glycanMass, x.Glycan), 1).tolist()


def load_from_file(path_to_file=mammalian_glycomedb_nlinked_path):
    gh = GlycanHypothesis(path_to_file)
    return gh.glycans


def pack_glycan_string(glycan_map, glycan_order):
    pack_string = "[{0}]".format(";".join([glycan_map[g] for g in glycan_order]))
    return pack_string


def parse_glycan_pack_string(pack_string, glycan_order):
    glycan_map = {g: int(c) for g, c in zip(glycan_order, re.findall('\\b\\d+\\b', pack_string))}
    return glycan_map
