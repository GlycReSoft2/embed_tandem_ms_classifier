from itertools import chain
from glycresoft_ms2_classification.ms import spectra
from glycresoft_ms2_classification.structure.sequence import Sequence


class namedfloat(float):
    '''Hack to shove a name into an already set __repr__ function'''
    def __new__(cls, value, name):
        obj = float.__new__(cls, value)
        obj.name = name
        return obj

    def __repr__(self):
        return "({} {})".format(self.name, self.real)

    def __str__(self):
        return repr(self)


def sequence_to_spectra(seq):
    if not isinstance(seq, Sequence):
        seq = Sequence(seq)
    precursor = spectra.ObservedPrecursorSpectrum([], [], 0, seq.mass, [])
    precursor.sequence = seq
    for fragment in chain.from_iterable(
            chain.from_iterable(
                zip(*map(seq.break_at, range(1, len(seq)))))):
        tandem = NamedTandemSpectrum(mass=fragment.get_mass(), charge=0,
                                     intensity=1, name=fragment.name)
        tandem.name = fragment.name
        precursor.tandem_data.append(tandem)
    precursor.tandem_data.sort(key=lambda x: x.neutral_mass)
    return precursor


class NamedTandemSpectrum(spectra.ObservedTandemSpectrum):
    def __repr__(self):
        rep = super(NamedTandemSpectrum, self).__repr__()
        return rep[:-1] + " {name}>".format(**self.other_data)
