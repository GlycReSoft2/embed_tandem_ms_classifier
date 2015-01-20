import re
from collections import Counter

import numpy as np
import pandas as pd

from glycresoft_ms2_classification.structure import sequence, constants
from glycresoft_ms2_classification.utils.collectiontools import groupby

fields_to_scan = [
    "bare_b_ions",
    "bare_y_ions",
    "b_ions_with_HexNAc",
    "y_ions_with_HexNAc"
]


key_patterns = {
    "bare_b_ions": re.compile(r"B(\d+)"),
    "bare_y_ions": re.compile(r"Y(\d+)"),
    "b_ions_with_HexNAc": re.compile(r"B(\d+)\+"),
    "y_ions_with_HexNAc": re.compile(r"Y(\d+)\+"),
}


def parse_key(key_string, ion_type):
    return int(key_patterns[ion_type].search(key_string).group(1)) - constants.FRAG_OFFSET


def tupelize(list_struct):
    return (list_struct[0], tuple(list_struct[1]),)


class BreakCounter(object):

    def __init__(self, source_sequence):
        try:
            self.sequence = sequence.Sequence(source_sequence)
        except:
            self.sequence = source_sequence
        self.ions_matched = set()
        self.b_breaks = pd.Series(np.zeros(len(self)))
        self.y_breaks = pd.Series(np.zeros(len(self)))
        self.b_hex_breaks = pd.Series(np.zeros(len(self)))
        self.y_hex_breaks = pd.Series(np.zeros(len(self)))

    def count_matches(self, matches):
        for match in matches:
            for ion in match["bare_b_ions"]:
                self.ions_matched.add(ion["key"])
                self.b_breaks[parse_key(ion["key"], "bare_b_ions")] += 1
            for ion in match["bare_y_ions"]:
                self.ions_matched.add(ion["key"])
                self.y_breaks[parse_key(ion["key"], "bare_y_ions")] += 1
            for ion in match["b_ions_with_HexNAc"]:
                self.ions_matched.add(ion["key"])
                self.b_hex_breaks[parse_key(ion["key"], "b_ions_with_HexNAc")] += 1
            for ion in match["y_ions_with_HexNAc"]:
                self.ions_matched.add(ion["key"])
                self.y_hex_breaks[parse_key(ion["key"], "y_ions_with_HexNAc")] += 1

    def estimate_break_parameter(self, lambda_range):
        pass

    def empty(self):
        num_matches = sum(self.b_breaks) + sum(self.y_breaks) + sum(self.b_hex_breaks) +\
            sum(self.y_hex_breaks)
        return not (num_matches > 0.0)

    def empty_hex(self):
        num_matches = sum(self.b_hex_breaks) + sum(self.y_hex_breaks)
        return not (num_matches > 0.0)

    def __repr__(self):
        return "b-ions: {0}\ny-ions: {1}\nb-hex: {2}\ny-hex: {3}\n{4}\n{5}".format(
            self.b_breaks.as_matrix(), self.y_breaks.as_matrix(),
            self.b_hex_breaks.as_matrix(), self.y_hex_breaks.as_matrix(), self.sequence, self.ions_matched)

    def __len__(self):
        if isinstance(self.sequence, sequence.Sequence):
            return len(self.sequence)
        else:
            return len(self.sequence[0])

    def __add__(self, other):
        if not isinstance(other, BreakCounter):
            raise TypeError("Cannot add a non-BreakCounter to a BreakCounter, {0}".format(type(other)))
        counter = BreakCounter((self.sequence, other.sequence))
        counter.b_breaks = self.b_breaks + other.b_breaks
        counter.y_breaks = self.y_breaks + other.y_breaks
        counter.b_hex_breaks = self.b_hex_breaks + other.b_hex_breaks
        counter.y_hex_breaks = self.y_hex_breaks + other.y_hex_breaks
        counter.ions_matched = self.ions_matched | other.ions_matched
        return counter

    def hist(self, kind="b", **kwargs):
        dim = None
        if kind == "b":
            dim = self.b_breaks
        elif kind == "y":
            dim = self.y_breaks
        elif kind == "bhex":
            dim = self.b_hex_breaks
        elif kind == "yhex":
            dim = self.y_hex_breaks

        return dim.plot(kind="bar", **kwargs)


def count_breaks(matched_ions):
    break_counts = dict()
    for glycopeptide, matches in groupby(matched_ions,
                                         lambda x: x["Glycopeptide_identifier"]).items():
        counter = BreakCounter(glycopeptide)
        counter.count_matches(matches)
        if counter.empty_hex():
            continue
        break_counts[glycopeptide] = counter
    return break_counts


class BondBreakageCounter(object):
    def __init__(self, *args, **kwargs):
        self.data = pd.Series()

    def __getitem__(self, key):
        try:
            self.data[key]
        except:
            return 0

    def __setitem__(self, key, value):
        try:
            self.data[key] = value
        except:
            pass

    def plot(self, *args, **kwargs):
        return self.data.plot(*args, **kwargs)

    def normalize(self):
        pass
