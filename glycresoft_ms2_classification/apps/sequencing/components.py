import logging
import itertools
from operator import attrgetter

from glycresoft_ms2_classification.structure.sequence import sequence_to_mass, Sequence, sequence_tokenizer
from glycresoft_ms2_classification.structure import residue, modification
from glycresoft_ms2_classification.structure.composition import Composition
from glycresoft_ms2_classification.ms import pickle

masser = attrgetter('mass')
precursor_mass_shift = Composition('H2O').mass
residue_symbols = dict(residue.symbol_to_residue)
residue_symbols.pop("I")

logging.basicConfig(level='DEBUG', format='%(asctime)s %(message)s', datefmt='%I:%M:%S %p')
logger = logging.getLogger()


precursor_mass_shift = Composition('H2O').mass
y_mass_shift = Composition('H2O').mass + Composition('H').mass - Composition('e').mass
b_mass_shift = Composition('H').mass - Composition('e').mass


def simple_mass(res, mods):
    mass = res.mass
    for mod in mods:
        mass += mod.mass
    return mass


class SequencePosition(tuple):
    def __new__(self, residue, modifications):
        mass = residue.mass
        mass += sum(m.mass for m in modifications)
        return tuple.__new__(self, (residue, tuple(modifications), mass))

    @property
    def mass(self):
        return self[2]

    def __len__(self):
        return 1


def generate_component_set(constant_modifications, variable_modifications):
    components = set(residue_symbols)
    const_modified = set()
    for mod in constant_modifications:
        for symbol, name in residue_symbols.items():
            if mod.valid_site(name, "internal"):
                components.add(symbol + "(" + mod.serialize() + ")")
                const_modified.add(symbol)

    # After applying constant modifications, we can't have the same residue unmodified, so
    # remove it from the component set.
    map(components.discard, const_modified)

    for mod in variable_modifications:
        for symbol, name in residue_symbols.items():
            if mod.valid_site(name):
                components.add(symbol + "(" + mod.serialize() + ")")

    components = list(components)
    return map(SimpleFragment, components)


def unordered_combinations(parts, n=2):
    for comb in itertools.combinations_with_replacement(parts, n):
        yield SimpleFragment(comb, mass=sum(masser(c) for c in comb), length=n, unordered=True)


class SqlSequenceBase(object):

    def to_sql(self):
        yield '''insert into Components
        (component_id, neutral_mass, serial) VALUES ({}, {}, {})'''.format(
            id, self.mass, pickle.dumps(self))


class SimpleFragment(SqlSequenceBase):
    def __init__(self, sequence, mass=None, missed_cleavages=0, length=None, unordered=False):
        if mass is None:
            if isinstance(sequence, basestring):
                mass = sequence_to_mass(str(sequence)) - precursor_mass_shift
                length = len(sequence_tokenizer(sequence)[0])
            else:
                mass = sum(map(masser, sequence))
                length = len(sequence)
        self.sequence = sequence
        self.mass = mass
        self.missed_cleavages = missed_cleavages
        self.length = length
        self.unordered = unordered

    def __iter__(self):
        return iter(self.sequence)

    def __len__(self):
        return self.length

    def __repr__(self):
        return "{ord}{seq}".format(ord="<U>" if self.unordered else "", seq=self.sequence)

    def to_sequence_positions(self):
        return [SequencePosition(r, m) for r, m in Sequence(self.sequence)]

    def to_sequence(self):
        if self.unordered:
            for perm in itertools.permutations(self.sequence):
                yield perm
        else:
            yield self.sequence


class SequenceCollection(SqlSequenceBase):
    def __init__(self, components):
        self.components = components
        self.mass = sum(c.mass for c in self.components)
        self.length = sum(len(c) for c in self.components)

    def __len__(self):
        return self.length

    def __repr__(self):
        return repr(self.components)

    def to_sequence(self):
        parts = [list(c.to_sequence()) for c in self.components]
        for path in itertools.product(*parts):
            yield ''.join(map(str, itertools.chain.from_iterable(path)))
