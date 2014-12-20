import os
import cPickle
import gzip
from collections import Iterable

from .utils import ParallelParser
from .spectra import (neutral_mass, matched_spectra,
                      ObservedPrecursorSpectra, ObservedTandemSpectra)


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
