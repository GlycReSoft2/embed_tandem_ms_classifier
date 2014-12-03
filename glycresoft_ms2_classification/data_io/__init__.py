import os
import cPickle
import gzip
from collections import Iterable

from ..structure.composition import composition_to_mass
PROTON = composition_to_mass("p")

MAX_CHARGE_STATE = 4
MACHINE_ACQUISITION_RANGE = 3200 * 24

opener = open
compressed_opener = gzip.open

opener_map = {
    "gz": compressed_opener,
    "": opener
}


def get_opener(file_name=""):
    ext = os.path.splitext(file_name)[1]
    return opener_map.get(ext, opener)


class DeconIOBase(object):

    @staticmethod
    def prepare(data_dict):
        holder = DeconIOBase(None)
        holder.data = data_dict
        return holder

    def __init__(self, file_path=None):
        if file_path is None:
            self.data = dict()
        else:
            self._load(file_path)

    def _load(self, file_path):
        self.data = cPickle.load(get_opener(file_path)(file_path, "rb"))

    def __iter__(self):
        return iter(self.data.items())

    def subset(self, indices):
        if not isinstance(indices, Iterable):
            indices = [indices]
        return self.prepare({i: self.data[i] for i in indices})

    def __getitem__(self, ix):
        return self.data[ix]


def default_loader(file_path):
    return cPickle.load(get_opener(file_path)(file_path, "rb"))


def default_serializer(obj, file_path):
    cPickle.dump(obj, get_opener(file_path)(file_path, "wb"))


def neutral_mass(mz, z):
    return (mz * z) - (z * PROTON)


class ObservedPrecursorSpectra(object):

    def __init__(self, scan_data, scan_ids, neutral_mass, tandem_data, **data):
        self.scans = scan_data
        self.scan_ids = scan_ids
        self.neutral_mass = neutral_mass
        self.tandem_data = tandem_data
        self.other_data = data
        self._iterkey = None

    def __repr__(self):
        return "<Scans: ({scans[0]}), Neutral Mass: {neutral_mass}, Tandem Spectra: {num_tandem}>".format(
            num_tandem=len(self.tandem_data),
            min_scan=self.scans[-1],
            **self.__dict__)

    def __ge__(self, o):
        return self.neutral_mass >= o

    def __gt__(self, o):
        return self.neutral_mass > o

    def __lt__(self, o):
        return self.neutral_mass < o

    def __le__(self, o):
        return self.neutral_mass <= o

    def __rge__(self, o):
        return self.neutral_mass <= o

    def __rgt__(self, o):
        return self.neutral_mass < o

    def __rlt__(self, o):
        return self.neutral_mass > o

    def __rle__(self, o):
        return self.neutral_mass >= o


class ObservedTandemSpectra(object):

    def __init__(self, mass, z, intensity, **data):
        self.charge = z
        self.mass = mass
        self.intensity = intensity
        self.other_data = data

    def __repr__(self):
        return "<ObservedTandemSpectra {neutral_mass}, {charge}, {intensity}>".format(**self.__dict__)

    def __ge__(self, o):
        return self.neutral_mass >= o

    def __gt__(self, o):
        return self.neutral_mass > o

    def __lt__(self, o):
        return self.neutral_mass < o

    def __le__(self, o):
        return self.neutral_mass <= o

    def __rge__(self, o):
        return self.neutral_mass <= o

    def __rgt__(self, o):
        return self.neutral_mass < o

    def __rlt__(self, o):
        return self.neutral_mass > o

    def __rle__(self, o):
        return self.neutral_mass >= o


def matched_spectra(scan_ids, neutral_mass, ppm_error, match_key, **data):
    return dict(scan_ids=scan_ids,
                neutral_mass=neutral_mass,
                ppm_error=ppm_error,
                key=match_key,
                other_data=data)
