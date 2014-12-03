import itertools
import pandas as pd


import classify_matches


def merge_predictions(*prediction_frames):
    for i, frame in enumerate(prediction_frames):
        frame["group_index"] = i
    merged_set = pd.concat(prediction_frames)
    all_matches = merged_set.groupby("Glycopeptide_identifier")
    for ident, matches in all_matches:
        pass
        # Merge common
        # Merge all


def do_merge_common(matches, n):
    if(len(matches.index) == n):
        pass

    return None


def do_combine_data(matches):
    pass


def group_ions_by_key(ion_series):
    pass
