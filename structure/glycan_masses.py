import csv
import os

from collections import OrderedDict
from operator import itemgetter

from structure.modification import AnonymousModificationRule

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

    def as_modification(self):
        return AnonymousModificationRule("Glycan" + self.composition, self.molecular_weight)

    molecular_weight = property(itemgetter(0), doc='Alias for field number 0')

    composition = property(itemgetter(1), doc='Alias for field number 1')


class GlycanHypothesis(list):
    def __init__(self, csv_path=None, records=None):
        if csv_path is not None:
            reader = csv.DictReader(open(csv_path))
            self.extend(reader)
        elif records is not None:
            self.extend(records)
        self.extract_glycans()

    def extract_glycans(self):
        self.glycans = list({SimpleGlycan(float(r["Molecular Weight"]), r["Compositions"]) for r in self})
        return self.glycans


def load_from_file(path_to_file=mammalian_glycomedb_nlinked_path):
    gh = GlycanHypothesis(path_to_file)
    return gh.glycans
