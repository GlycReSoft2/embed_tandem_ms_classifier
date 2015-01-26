from collections import Counter
from copy import deepcopy

import sqlitedict
import pandas as pd
import numpy as np

from ..data_processing import prepare_model_file
from ...structure import sequence
from ...ms import default_loader, spectra


def take_best_predictions_spectra(predictions, spectra_db):
    threshold_predictions = predictions.kvquery()
    collector = sqlitedict.open("./spectra_counting.db", tablename="spectra_map")
    for ix, pred in threshold_predictions.iterrows():
        glycopeptide_ident = pred.Glycopeptide_identifier
        print(glycopeptide_ident)
        matched_spectra = [spectra_db[scan_id] for scan_id in pred.scan_id_range]
        collector[ix] = {"glycopeptide": glycopeptide_ident, "spectra": matched_spectra}
    collector.commit()
    return collector


def sequence_position_to_tuple(seq_pos):
    mods = tuple(seq_pos[1])
    pos = tuple(seq_pos[0], mods)
    return pos


def bond_pair(bond_pair):
    return frozenset(map(sequence_position_to_tuple, bond_pair))


def interpolate_key_bond_map(sequence_obj):
    key_bond_map = {}
    for i in range(1, len(sequence_obj)):
        b_ions, y_ions = sequence_obj.break_at(i)
        bond = bond_pair(sequence_obj[(i-1):(i+1)])
        for f in b_ions:
            key_bond_map[f.name] = bond
        for f in y_ions:
            key_bond_map[f.name] = bond
    return key_bond_map


def main(predictions, matched_spectra_db):
    matched_spectra = take_best_predictions_spectra(predictions, matched_spectra_db)
    counter = Counter()
    for glycopeptide, spectra in matched_spectra.iteritems():
        bond_map = interpolate_key_bond_map(sequence.Sequence(glycopeptide))
        for spectrum in spectra:
            for tandem in spectrum.tandem_data:
                print(len(tandem.annotation))