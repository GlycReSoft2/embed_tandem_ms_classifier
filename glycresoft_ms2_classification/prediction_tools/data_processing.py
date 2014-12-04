import json
from math import fabs
import re
from cStringIO import StringIO

import pandas as pd
import numpy as np

from .sequence_handling import get_sequence_length
from .sequence_handling import modification_signature
from .sequence_handling import glycosites
from .sequence_handling import percent_expected_ions_with_hexnac_observed

from ..structure import sequence
from ..structure import modification

# Defines the list of fields that must be JSON encoded and decoded
# for storage in CSV files. Because regular scientists cannot possibly
# understand the fractal madness of pure hierarchical structured data
# files. YAML is out of the question.
json_serialized_columns = ["Oxonium_ions", "Stub_ions",
                           "b_ion_coverage", "b_ions_with_HexNAc",
                           "y_ion_coverage", "y_ions_with_HexNAc",
                           "peptideCoverageMap", "hexNAcCoverageMap",
                           "scan_id_range"]


# Decode :json_serialized_columns columns into Python objects for interaction.
# These columns are not "fast", and should only be accessed for
# big tasks. This function destructively modifies the input object.
def deserialize_compound_fields(data_struct):
    for field in json_serialized_columns:
        try:
            data_struct[field] = pd.Series(map(json.loads, data_struct[field]))
        except TypeError:
            #print("Deserializing %s failed (TypeError)" % field)
            pass
        except KeyError:
            #print("Deserializing %s failed (KeyError)" % field)
            pass
        except ValueError, e:
            raise ValueError(str(e) + " " + field)

    return data_struct


# Encode :json_serialized_columns columns into JSON strings for writing
# to CSV. This is a destructive operation and should be called on a copy
# of the working object, or else the object will need to go through
# :deserialize_compound_fields again. This function destructively modifies
# the input object.
def serialize_compound_fields(data_struct):
    for field in json_serialized_columns:
        try:
            data_struct[field] = map(json.dumps, data_struct[field])
        except:
            pass  # Not all fields must be present
    return data_struct


# If the file has not had its percentage coverage by ion family statistics
# computed, but has all of the necessary information, compute all the ion family
# statistics and store them destructively. This function destructively
# modifies the input object.
def shim_percent_calcs(data):
    if "percent_b_ion_coverage" not in data:
        if "total_b_ions_possible" not in data:
            data["percent_b_ion_coverage"] = 0
        else:
            data["percent_b_ion_coverage"] = map(
                len, data["b_ion_coverage"]) / data['total_b_ions']
    if "percent_y_ion_coverage" not in data:
        if "total_y_ions_possible" not in data:
            data["percent_y_ion_coverage"] = 0
        else:
            data["percent_y_ion_coverage"] = map(
                len, data["y_ion_coverage"]) / data['total_y_ions']

    if "percent_b_ion_with_HexNAc_coverage" not in data:
        if "possible_b_ions_HexNAc" not in data:
            data["percent_b_ion_with_HexNAc_coverage"] = 0
        else:
            data["percent_b_ion_with_HexNAc_coverage"] = map(
                len, data["possible_b_ions_HexNAc"]) / data['possible_b_ions_HexNAc']

    if "percent_y_ion_with_HexNAc_coverage" not in data:
        if "possible_y_ions_HexNAc" not in data:
            data["percent_y_ion_with_HexNAc_coverage"] = 0
        else:
            data["percent_y_ion_with_HexNAc_coverage"] = map(
                len, data["possible_y_ions_HexNAc"]) / data['possible_y_ions_HexNAc']


# Given the path to a CSV file for a model or target to be classified, parse the
# CSV into a pandas.DataFrame and make sure it has all of the necessary columns and
# transformations to either build a model with it or
def prepare_model_file(path):
    if path[-4:] == "json":
        return PredictionResults.deserialize(path)
    data = pd.read_csv(path)
    deserialize_compound_fields(data)
    data.Peptide_mod = data.Peptide_mod.astype(str).replace("nan", "")
    data.ix[np.isnan(data.vol)].vol = -1
    data.ix[np.isnan(data.startAA)].startAA = -1
    data.ix[np.isnan(data.endAA)].endAA = -1
    data['abs_ppm_error'] = data.ppm_error.abs()
    if "peptideLens" not in data:
        data = get_sequence_length(data)
    if "numOxIons" not in data:
        data["numOxIons"] = map(len, data["Oxonium_ions"])
    if "numStubs" not in data:
        data["numStubs"] = map(len, data["Stub_ions"])
    if "percent_b_ion_coverage" not in data:
        shim_percent_calcs(data)
    data["glycosylation_sites"] = data.Glycopeptide_identifier.apply(glycosites)
    ion_coverage_score(data)
    modification_signature(data)
    infer_glycan_mass(data)

    return data


# Serialize the model object into a CSV file. Creates a copy of the data which has its
# complex fields JSON encoded.
def save_model_file(data_struct, path):
    if isinstance(data_struct, PredictionResults):
        data_struct.serialize(path)
        return
    write_copy = data_struct.copy()
    serialize_compound_fields(write_copy)
    write_copy.to_csv(path, index=False, encoding="utf-8")


def ion_coverage_score(data_struct):
    coverage_data = data_struct.apply(compute_ion_coverage_map, 1)
    data_struct['meanCoverage'] = coverage_data["meanCoverage"]
    data_struct['percentUncovered'] = coverage_data["percentUncovered"]
    data_struct['meanHexNAcCoverage'] = coverage_data['meanHexNAcCoverage']
    data_struct['peptideCoverageMap'] = coverage_data['peptideCoverageMap']
    data_struct['hexNAcCoverageMap'] = coverage_data['hexNAcCoverageMap']
    data_struct["bIonCoverageMap"] = coverage_data["bIonCoverage"]
    data_struct["bIonCoverageWithHexNAcMap"] = coverage_data[
        "bIonCoverageWithHexNAc"]
    data_struct["yIonCoverageMap"] = coverage_data["yIonCoverage"]
    data_struct["yIonCoverageWithHexNAcMap"] = coverage_data[
        "yIonCoverageWithHexNAc"]
    return data_struct


# Maps the coverage along the peptide backbone by matching ions
# of the correct size to the backbone.
def compute_ion_coverage_map(row):
    peptide_length = row['peptideLens']
    total_coverage = np.zeros(peptide_length)
    # Could be made more efficient by building by matrix operations

    b_ion_coverage = np.zeros(peptide_length)
    for b_ion in row['b_ion_coverage']:
        ion = np.zeros(peptide_length)
        ion[:int(b_ion['key'].replace("B", ''))] = 1
        b_ion_coverage += ion

    y_ion_coverage = np.zeros(peptide_length)
    for y_ion in row['y_ion_coverage']:
        ion = np.zeros(peptide_length)
        ion[:int(y_ion['key'].replace("Y", ''))] = 1
        y_ion_coverage += ion

    b_ions_with_HexNAc = np.zeros(peptide_length)
    for b_ion in row['b_ions_with_HexNAc']:
        ion = np.zeros(peptide_length)
        ion[:int(re.findall(r'B(\d+)\+', b_ion['key'])[0])] = 1
        b_ions_with_HexNAc += ion

    y_ions_with_HexNAc = np.zeros(peptide_length)
    for y_ion in row['y_ions_with_HexNAc']:
        ion = np.zeros(peptide_length)
        ion[:int(re.findall(r'Y(\d+)\+', y_ion['key'])[0])] = 1
        y_ions_with_HexNAc += ion

    total_coverage += b_ion_coverage.astype(
        np.int32) | b_ions_with_HexNAc.astype(np.int32)
    total_coverage += y_ion_coverage[
        ::-1].astype(np.int32) | y_ions_with_HexNAc[::-1].astype(np.int32)  # Reverse

    b_ions_covered = sum((b_ion_coverage.astype(
        np.int32) | b_ions_with_HexNAc.astype(np.int32)) != 0)

    y_ions_covered = sum((y_ion_coverage[
        ::-1].astype(np.int32) | y_ions_with_HexNAc[::-1].astype(np.int32)) != 0)

    expected_ions = (peptide_length * 2.0) - 1
    percent_expected_observed = (b_ions_covered + y_ions_covered) / expected_ions

    percent_uncovered = (sum(total_coverage == 0) / float(peptide_length))
    mean_coverage = total_coverage.mean() / (float(peptide_length))

    hexnac_coverage = (b_ions_with_HexNAc + y_ions_with_HexNAc[::-1])/2.

    mean_hexnac_coverage = percent_expected_ions_with_hexnac_observed(row)

    return pd.Series({"meanCoverage": mean_coverage,
                      "percentUncovered": percent_uncovered,
                      "percentExpectedObserved": percent_expected_observed,
                      "meanHexNAcCoverage": mean_hexnac_coverage,
                      "peptideCoverageMap": list(total_coverage),
                      "hexNAcCoverageMap": list(hexnac_coverage),
                      "bIonCoverage": b_ion_coverage + b_ions_with_HexNAc,
                      "bIonCoverageWithHexNAc": b_ions_with_HexNAc,
                      "yIonCoverage": y_ion_coverage[::-1] + y_ions_with_HexNAc[::-1],
                      "yIonCoverageWithHexNAc":  y_ions_with_HexNAc[::-1]
                      })


def build_ion_map(ion_set, ion_type, length):
    pat = re.compile(r"{0}(\d+)\+?".format(ion_type))
    coverage = np.zeros(length)
    for ion in ion_set:
        inst_cover = np.zeros(length)
        inst_cover[:int(pat.findall(ion["key"])[0])] = 1
        coverage += inst_cover

    return coverage


# Based upon the assumption that the same MS1 Score + Observed Mass implies
# the same precursor ion, the data is grouped by these two scores, and sort within groups.
# If the input DataFrame does not have an MS2 Score column, it will be removed
# from the within-group sorting criteria.
def determine_ambiguity(data_struct):
    data_struct['ambiguity'] = False
    groupings = data_struct.groupby(["MS1_Score", "Obs_Mass"])
    # Bit list for pandas.DataFrame sort function
    data_struct["_abs_ppm_error"] = data_struct["ppm_error"].abs()
    sort_direction = [1, 1, 1, 0, 1, 0, 1]
    sort_fields = ["MS2_Score", "numStubs", "meanCoverage",
                   "percentUncovered", "peptideLens",
                   "_abs_ppm_error", "meanHexNAcCoverage"]
    if "MS2_Score" not in data_struct:
        sort_direction.pop(0)
        sort_fields.pop(0)
    regrouping = []
    for name, group in groupings:
        if (len(group.index) == 1):
            regrouping.append(group)
            continue
        group = group.sort(columns=sort_fields, ascending=sort_direction,
                           axis=0)
        group['ambiguity'] = True
        regrouping.append(group)

    regroup = pd.concat(regrouping)
    regroup.drop("_abs_ppm_error", axis=1, inplace=True)
    return regroup


# A comparison function based upon the sorting method used in
# :determine_ambiguity
def compare_rows(a, b):
    sort_direction = [1, 1, 1, 0, 1, 0, 1]
    sort_fields = ["MS2_Score", "numStubs", "meanCoverage",
                   "percentUncovered", "peptideLens",
                   "abs_ppm_error", "meanHexNAcCoverage"]
    a_score = 0
    b_score = 0
    for field_i, field in enumerate(sort_fields):
        if a[field] > b[field]:
            if sort_direction[field_i]:
                a_score += 1
            else:
                b_score += 1
    return a_score - b_score


def call_by_coverage(data_struct, criterion_fn=None):
    call = None
    if criterion_fn is None:
        call = np.logical_and(data_struct.meanHexNAcCoverage >
                              0.2, data_struct.percentUncovered < 0.2)
    else:
        call = data_struct.apply(criterion_fn, 1)
    return call


def infer_glycan_mass(data_struct):
    hexNAc_mass = modification.ModificationTable.bootstrap()["HexNAc"].mass
    glycopeptide_mass = data_struct.Glycopeptide_identifier.apply(lambda x: sequence.Sequence(x).mass)
    glycan_mass = data_struct.Calc_mass - (glycopeptide_mass - (data_struct.glyco_sites * hexNAc_mass))
    data_struct["glycanMass"] = glycan_mass
    return data_struct


#TODO
def recalibrate_mass_accuracy(data_struct, ms1_tolerance, ms2_tolerance):
    ms1_recalibrated = data_struct.ppm_error.abs() <= ms1_tolerance
    recalibrated = data_struct.ix[ms1_recalibrated].apply(
        filter_tandem_matches_by_ppm, 1, args=(ms2_tolerance))
    ion_coverage_score(recalibrated)
    recalibrated.noise_filter = None
    recalibrated.MS2_Score = None


def filter_tandem_matches_by_ppm(row, ms2_tolerance):
    for ion_type in ["Oxonium_ions", "Stub_ions",
                     "b_ion_coverage", "b_ions_with_HexNAc",
                     "y_ion_coverage", "y_ions_with_HexNAc"]:
        row[ion_type] = [
            ion for ion in row[ion_type] if fabs(ion['ppm_error']) <= ms2_tolerance]
    return row


class PredictionResultsJSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, pd.DataFrame):
            return json.loads(o.to_json(orient="records"))
        else:
            return json.JSONEncoder.default(self, o)


class PredictionResults(object):
    def __init__(self, metadata, predictions):
        self.metadata = metadata
        self.predictions = predictions

    def __getitem__(self, inds):
        return self.predictions[inds]

    def __setitem__(self, inds, value):
        self.predictions[inds] = value

    def __getattr__(self, name):
        try:
            return object.__getattribute__(self, name)
        except:
            try:
                return getattr(object.__getattribute__(self, "predictions"), name)
            except:
                raise AttributeError("Prediction Results has no attribute {attr}".format(attr=name))

    def __contains__(self, name):
        return name in self.predictions

    @classmethod
    def deserialize(cls, data_buffer):
        if isinstance(data_buffer, basestring):
            data_buffer = open(data_buffer, "rb")
        loose_data = json.load(data_buffer)
        metadata = loose_data["metadata"]
        predictions = loose_data["predictions"]
        predictions = cls.ensure_fields(predictions)
        instance = cls(metadata, predictions)
        return instance

    @classmethod
    def prepare(cls, metadata, predictions):
        prediction_frame = PredictionResults(metadata, predictions)
        prediction_frame.predictions = prediction_frame.ensure_fields(predictions)
        return prediction_frame

    def serialize(self, target_buffer=None):
        if target_buffer is None:
            target_buffer = StringIO()
        if isinstance(target_buffer, basestring):
            target_buffer = open(target_buffer, "wb")
        encoder = PredictionResultsJSONEncoder(encoding="ascii")
        for chunk in encoder.iterencode(self.__dict__):
            target_buffer.write(chunk)
        return target_buffer

    @staticmethod
    def ensure_fields(raw_predictions):
        data = pd.DataFrame(raw_predictions)
        data.Peptide_mod = data.Peptide_mod.astype(str).replace("nan", "")
        # Handle Decoys
        data.ix[np.isnan(data.vol)].vol = -1
        data.ix[np.isnan(data.startAA)].startAA = -1
        data.ix[np.isnan(data.endAA)].endAA = -1
        if "bad_oxonium_ions" not in data:
            data["bad_oxonium_ions"] = data.Oxonium_ions.apply(lambda x: [])
        if "peptideLens" not in data:
            data = get_sequence_length(data)
        if "numOxIons" not in data:
            data["numOxIons"] = map(len, data["Oxonium_ions"])
        if "numStubs" not in data:
            data["numStubs"] = map(len, data["Stub_ions"])
        if "percent_b_ion_coverage" not in data:
            shim_percent_calcs(data)
        data["glycosylation_sites"] = data.Glycopeptide_identifier.apply(glycosites)
        ion_coverage_score(data)
        modification_signature(data)
        infer_glycan_mass(data)
        return data


def convert_csv_to_nested(path):
    frame = prepare_model_file(path)
    metadata = {
        "ms1_results_file": path,
        "_converted": True
    }
    pr = PredictionResults(metadata, frame)
    return pr
