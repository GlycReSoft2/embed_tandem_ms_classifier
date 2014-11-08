import itertools
import classify_matches

import pandas as pd


def merge_predictions(*prediction_sets):
    merged_set = pd.concat(prediction_sets)
    all_matches = merged_set.groupby("Glycopeptide_identifier")
    common_matches = []
    merged_matches = []
    n_sets = len(prediction_sets)
    for ident, matches in all_matches:
        pass
        # Merge common
        # Merge all

def do_merge_common(matches, n):
    if(len(group.index) == n):
        pass

    return None

def do_combine_data(matches):
    pass
