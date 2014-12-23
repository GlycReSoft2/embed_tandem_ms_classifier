import itertools
from math import fabs

import pandas as pd

import prediction_tools
from . import collectiontools


def merge_ion_matches(matches):
    groups = collectiontools.groupby(itertools.chain.from_iterable(matches),
                                     lambda x: x["key"])
    best_matches = []
    for key, matched_key in groups.items():
        best_match = matched_key[0]
        best_ppm = fabs(best_match["ppm_error"])
        ppm_to_scan_id = {best_match["scan_id"]: best_match["ppm_error"]}
        for match in matched_key[1:]:
            ppm_to_scan_id[match["scan_id"]] = match["ppm_error"]
            if fabs(match["ppm_error"]) < best_ppm:
                best_match = match
                best_ppm = fabs(match["ppm_error"])
        best_match["scan_map"] = ppm_to_scan_id
        best_matches.append(best_match)
    return best_matches


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
