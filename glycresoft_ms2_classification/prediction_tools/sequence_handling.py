import re

from collections import defaultdict

import numpy as np
import pandas as pd

from ..structure.sequence import sequence_tokenizer, sequence_length, Sequence, golden_pair_map
from ..structure import stub_glycopeptides
from ..structure import constants as structure_constants

StubGlycopeptide = stub_glycopeptides.StubGlycopeptide


def get_sequence_length(data):
    data["peptideLens"] = data.Peptide.apply(sequence_length)
    return data


def extract_modifications(sequence):
    sequence_tokens, mods, glycan, n_term, c_term = sequence_tokenizer(sequence)
    mods = sorted(mods)
    return ','.join(mods) + glycan


def modification_signature(data_struct):
    mod_data = data_struct.Glycopeptide_identifier.apply(extract_modifications)
    data_struct["modificationSignature"] = mod_data
    return data_struct


def find_occupied_glycosylation_sites(sequence):
    tokens, mods, glycan, n_term, c_term = sequence_tokenizer(sequence)
    glycosites = [i for i in range(len(tokens)) if "HexNAc" in tokens[i][1]]
    return glycosites


def total_expected_ions_with_hexnac(row, kind):
    if kind not in ["b", "y"]:
        raise KeyError("Can calculate total_expected_ions_with_hexnac\
         for b or y ions. {kind} requested".format(kind=kind))
    if kind == "b":
        # exclude b1 ion
        return row.peptideLens - min(row.glycosylation_sites) - structure_constants.EXCLUDE_B1
    else:
        # Distance from the end to last glycosylation site
        return row.peptideLens - (row.peptideLens - max(row.glycosylation_sites))


def percent_expected_ions_with_hexnac_observed(row):
    expected_b_hexnac = total_expected_ions_with_hexnac(row, "b")
    expected_y_hexnac = total_expected_ions_with_hexnac(row, "y")
    b_ions_hexnac_observed = len(row.b_ions_with_HexNAc)
    y_ions_hexnac_observed = len(row.y_ions_with_HexNAc)

    percent_expected_observed = (b_ions_hexnac_observed + y_ions_hexnac_observed) \
        / float(expected_b_hexnac + expected_y_hexnac)

    return percent_expected_observed


# Maps the coverage along the peptide backbone by matching ions
# of the correct size to the backbone.
def compute_ion_coverage_map(glycopeptide_match, include_golden_pairs=True):
    """Compute the per-ion type coverage for the current glycopeptide
    and integrate complementar ion ladders to compute the mean per residue
    coverage with and without HexNAc.

    Parameters
    ----------
    glycopeptide_match: pandas.Series
        Observed ion mapping for a single glycopeptide prediction
    include_golden_pairs: bool
        Also compute observed golden pair ions

    Returns
    -------
    pandas.Series:
        Coverage statistics and per site coverage maps
    """
    peptide_length = glycopeptide_match['peptideLens']
    total_coverage = np.zeros(peptide_length)

    b_ion_coverage = np.zeros(peptide_length)
    for b_ion in glycopeptide_match['b_ion_coverage']:
        ion = np.zeros(peptide_length)
        ion[:int(b_ion['key'].replace("B", '')) - 1] = 1
        b_ion_coverage += ion

    y_ion_coverage = np.zeros(peptide_length)
    for y_ion in glycopeptide_match['y_ion_coverage']:
        ion = np.zeros(peptide_length)
        ion[:int(y_ion['key'].replace("Y", '')) - 1] = 1
        y_ion_coverage += ion

    b_ions_with_HexNAc = np.zeros(peptide_length)
    for b_ion in glycopeptide_match['b_ions_with_HexNAc']:
        ion = np.zeros(peptide_length)
        ion[:int(re.findall(r'B(\d+)\+', b_ion['key'])[0]) - 1] = 1
        b_ions_with_HexNAc += ion

    y_ions_with_HexNAc = np.zeros(peptide_length)
    for y_ion in glycopeptide_match['y_ions_with_HexNAc']:
        ion = np.zeros(peptide_length)
        ion[:int(re.findall(r'Y(\d+)\+', y_ion['key'])[0]) - 1] = 1
        y_ions_with_HexNAc += ion

    # Combine the b and y ion series to compute the total coverage
    total_coverage += b_ion_coverage.astype(
        np.int32) | b_ions_with_HexNAc.astype(np.int32)
    total_coverage += y_ion_coverage[
        ::-1].astype(np.int32) | y_ions_with_HexNAc[::-1].astype(np.int32)  # Reverse

    b_ions_covered = sum((b_ion_coverage.astype(
        np.int32) | b_ions_with_HexNAc.astype(np.int32)) != 0)

    y_ions_covered = sum((y_ion_coverage[
        ::-1].astype(np.int32) | y_ions_with_HexNAc[::-1].astype(np.int32)) != 0)

    expected_ions = (peptide_length * 2.0) - 1
    percent_expected_observed = (
        b_ions_covered + y_ions_covered) / expected_ions

    percent_uncovered = (sum(total_coverage == 0) / float(peptide_length))
    mean_coverage = total_coverage.mean() / (float(peptide_length))

    hexnac_coverage = (b_ions_with_HexNAc + y_ions_with_HexNAc[::-1]) / 2.

    mean_hexnac_coverage = percent_expected_ions_with_hexnac_observed(glycopeptide_match)

    fragment_golden_pair_map = golden_pair_map(glycopeptide_match.Glycopeptide_identifier)
    fragments_observed = {f["key"]: f for f in glycopeptide_match['b_ions_with_HexNAc'] +
                          glycopeptide_match['y_ions_with_HexNAc'] + glycopeptide_match['b_ion_coverage'] +
                          glycopeptide_match['y_ion_coverage']}

    if include_golden_pairs:
        # Find each pair of fragments that would result from a single fragmentation
        golden_pairs_observed = 0
        golden_pairs_expected = 0
        golden_pairs_hexnac_observed = 0
        golden_pairs_hexnac_expected = 0
        [frag.pop("observed_golden_pairs", False) for frag in fragments_observed.values()]

        for key, frag in fragments_observed.items():
            golden_pairs = fragment_golden_pair_map[key]
            frag['golden_pairs'] = golden_pairs
            if "observed_golden_pairs" not in frag:
                frag["observed_golden_pairs"] = set()
            for golden_pair in golden_pairs:
                if golden_pair in frag["observed_golden_pairs"]:
                    continue
                if golden_pair in fragments_observed:
                    golden_pairs_observed += 1
                    if "HexNAc" in key or "HexNAc" in golden_pair:
                        golden_pairs_hexnac_observed += 1
                    frag["observed_golden_pairs"].add(golden_pair)
                    golden_pair_frag = fragments_observed[golden_pair]
                    if "observed_golden_pairs" not in golden_pair_frag:
                        golden_pair_frag["observed_golden_pairs"] = set()
                    golden_pair_frag["observed_golden_pairs"].add(key)
                golden_pairs_expected += 1
                if "HexNAc" in golden_pair or "HexNAc" in key:
                    golden_pairs_hexnac_expected += 1

        # Convert each set into a list
        for frag in fragments_observed.values():
            if "observed_golden_pairs" in frag:
                frag['observed_golden_pairs'] = list(frag['observed_golden_pairs'])

        golden_pairs_results = {
            "golden_pairs_observed": golden_pairs_observed,
            "golden_pairs_hexnac_observed": golden_pairs_hexnac_observed,
            "golden_pairs_expected": golden_pairs_expected,
            "golden_pairs_hexnac_expected": golden_pairs_hexnac_expected,
        }

    result = pd.Series({
        "meanCoverage": mean_coverage,
        "percentUncovered": percent_uncovered,
        "percentExpectedObserved": percent_expected_observed,
        "meanHexNAcCoverage": mean_hexnac_coverage,
        "peptideCoverageMap": list(total_coverage),
        "hexNAcCoverageMap": list(hexnac_coverage),
        "bIonCoverage": b_ion_coverage,
        "bIonCoverageWithHexNAc": b_ions_with_HexNAc,
        "yIonCoverage": y_ion_coverage[::-1],
        "yIonCoverageWithHexNAc":  y_ions_with_HexNAc[::-1],
    })

    if include_golden_pairs:
        result.update(golden_pairs_results)

    return result


def build_ion_map(ion_set, ion_type, length):
    pat = re.compile(r"{0}(\d+)\+?".format(ion_type))
    coverage = np.zeros(length)
    for ion in ion_set:
        inst_cover = np.zeros(length)
        inst_cover[:int(pat.findall(ion["key"])[0])] = 1
        coverage += inst_cover

    return coverage


def stubs_observed_expected_ratio(row):
    obs = row.numStubs
    expected = len(StubGlycopeptide.from_sequence(
        Sequence(row.Glycopeptide_identifier)).get_stubs())
    return obs/float(expected)
