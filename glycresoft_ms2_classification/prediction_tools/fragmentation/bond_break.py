import re
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
          matched_spectra = [spectra_db[scan_id] for scan_id in pred.scan_id_range]
          collector[ix] = {"glycopeptide": glycopeptide_ident, "spectra": matched_spectra}
      collector.commit()
      return collector

def sequence_position_to_tuple(seq_pos):
    mods = tuple(seq_pos[1])
    pos = tuple((seq_pos[0], mods))
    return pos


def bond_pair(bond_pair):
    return frozenset( ( tuple(map(sequence_position_to_tuple, bond_pair)  ), ) )

def interpolate_key_bond_map(sequence_obj):
    if not isinstance(sequence_obj, sequence.Sequence):
        sequence_obj = sequence.Sequence(sequence_obj)
    key_bond_map = {}
    for i in range(1, len(sequence_obj)):
        b_ions, y_ions = sequence_obj.break_at(i)
        bond = bond_pair(sequence_obj[(i-1):(i+1)])
        if len(bond) != 1:
            raise Exception(str(sequence_obj), bond, i-1, i+1)
        for f in b_ions:
            key_bond_map[f.name] = bond
        for f in y_ions:
            key_bond_map[f.name] = bond
    return key_bond_map


def fit(predictions, matched_spectra_db):
    matched_spectra = take_best_predictions_spectra(predictions, matched_spectra_db)
    counter = Counter()
    for psm in matched_spectra.itervalues():
        glycopeptide = psm["glycopeptide"]
        matched = psm["spectra"]
        bond_map = interpolate_key_bond_map(sequence.Sequence(glycopeptide))
        for spectrum in matched:
            for tandem in spectrum.tandem_data:
                if len(tandem.annotation) > 0:
                    for annot in tandem.annotation:
                        if annot[0] == glycopeptide and is_backbone_ion(annot[1]):
                            if len(bond_map[annot[1]]) != 1:
                                raise Exception(annot, bond_map[annot[1]], spectrum.scan_ids)
                            counter[bond_map[annot[1]]] += 1

    counter = pd.Series(counter)
    log_probabilities = np.log(counter / counter.sum())
    return log_probabilities

def is_backbone_ion(key):
    return re.match(r"^[BYby]\d+", key)

def evaluate(model, spectrum, seq):
    model_score = 0.0
    null_model_score = 0.0
    bond_map = interpolate_key_bond_map(seq)
    null_prob = np.log(1.0/len(model))
    cntr = 0
    for tandem_spectrum in spectrum.tandem_data:
        for case in tandem_spectrum.annotation:
            if case[0] == seq:
                ion_type = case[1]
                if ion_type in bond_map:
                    bond = bond_map[ion_type]
                    if bond in model:
                        model_score += model[bond]
                    else:
                        model_score += null_prob
                    cntr += 1
                    null_model_score += null_prob
    diff = model_score - null_model_score
    diff *= (cntr/float(len(bond_map)))
    diff = np.exp(diff)
    return (diff)
