import json
from math import fabs
from cStringIO import StringIO
from types import MethodType
from functools import partial

import sqlitedict
import pandas as pd
import numpy as np

from .sequence_handling import get_sequence_length, compute_ion_coverage_map
from .sequence_handling import modification_signature
from .sequence_handling import find_occupied_glycosylation_sites
from .sequence_handling import stubs_observed_expected_ratio

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
            print("Deserializing %s failed (TypeError)" % field)
            pass
        except KeyError:
            print("Deserializing %s failed (KeyError)" % field)
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
def prepare_model_file(path, recompute=False):
    if path[-4:] == "json":
        return PredictionResults.deserialize(path, recompute=recompute)
    data = pd.read_csv(path)
    deserialize_compound_fields(data)
    data.Peptide_mod = data.Peptide_mod.astype(str).replace("nan", "")
    data.ix[np.isnan(data.vol)].vol = -1
    data.ix[np.isnan(data.startAA)].startAA = -1
    data.ix[np.isnan(data.endAA)].endAA = -1
    data['abs_ppm_error'] = data.ppm_error.abs()
    if "peptideLens" not in data or recompute:
        data = get_sequence_length(data)
    if "numOxIons" not in data or recompute:
        data["numOxIons"] = map(len, data["Oxonium_ions"])
    if "numStubs" not in data or recompute:
        data["numStubs"] = map(len, data["Stub_ions"])
    if "percent_b_ion_coverage" not in data or recompute:
        shim_percent_calcs(data)
    if "glycosylation_sites" not in data or recompute:
        data["glycosylation_sites"] = data.Glycopeptide_identifier.apply(find_occupied_glycosylation_sites)
    if "meanCoverage" not in data or recompute:
        ion_coverage_score(data)
    if "modificationSignature" not in data or recompute:
        modification_signature(data)
    if "glycanMass" not in data or recompute:
        infer_glycan_mass(data)

    return data


# Serialize the model object into a CSV file. Creates a copy of the data which has its
# complex fields JSON encoded.
def save_model_file(data_struct, path):
    if isinstance(data_struct, PredictionResults):
        data_struct.serialize(path + ".json")
        return path + ".json"
    write_copy = data_struct.copy()
    serialize_compound_fields(write_copy)
    write_copy.to_csv(path + ".csv", index=False, encoding="utf-8")
    return path + ".csv"


def ion_coverage_score(data_struct):
    coverage_data = data_struct.apply(compute_ion_coverage_map, 1)
    data_struct['meanCoverage'] = coverage_data["meanCoverage"]
    data_struct['percentUncovered'] = coverage_data["percentUncovered"]
    data_struct['meanHexNAcCoverage'] = coverage_data['meanHexNAcCoverage']
    data_struct['peptideCoverageMap'] = coverage_data['peptideCoverageMap']
    data_struct['hexNAcCoverageMap'] = coverage_data['hexNAcCoverageMap']
    data_struct["bIonCoverageMap"] = coverage_data["bIonCoverage"]
    data_struct["bIonCoverageWithHexNAcMap"] = coverage_data["bIonCoverageWithHexNAc"]
    data_struct["yIonCoverageMap"] = coverage_data["yIonCoverage"]
    data_struct["yIonCoverageWithHexNAcMap"] = coverage_data["yIonCoverageWithHexNAc"]
    data_struct["goldenPairsObserved"] = coverage_data["golden_pairs_observed"]
    data_struct["goldenPairsHexNAcObserved"] = coverage_data["golden_pairs_hexnac_observed"]
    data_struct["goldenPairsExpected"] = coverage_data["golden_pairs_expected"]
    data_struct["goldenPairsHexNAcExpected"] = coverage_data["golden_pairs_hexnac_expected"]
    return data_struct


# Based upon the assumption that the same MS1 Score + Observed Mass implies
# the same precursor ion, the data is grouped by these two scores, and sort within groups.
# If the input DataFrame does not have an MS2 Score column, it will be removed
# from the within-group sorting criteria.
def determine_ambiguity(data_struct):
    wrap = None
    if isinstance(data_struct, PredictionResults):
        wrap = data_struct
        data_struct = data_struct.predictions
    data_struct['ambiguity'] = False
    groupings = data_struct.groupby(["MS1_Score", "Obs_Mass"])
    # Bit list for pandas.DataFrame sort function
    sort_direction = [1, 1, 1, 0, 1, 0, 1]
    sort_fields = ["MS2_Score", "numStubs", "meanCoverage",
                   "percentUncovered", "peptideLens",
                   "ppm_error", "meanHexNAcCoverage"]
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
    if wrap is not None:
        wrap.predictions = regroup
        return wrap
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
    glycopeptide_mass = data_struct.Glycopeptide_identifier.apply(
        lambda x: sequence.Sequence(x).mass)
    glycan_mass = data_struct.Calc_mass - \
        (glycopeptide_mass - (data_struct.glyco_sites * hexNAc_mass))
    data_struct["glycanMass"] = glycan_mass
    return data_struct


# TODO
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


def query_threshold(MS2_Score=0, meanCoverage=0, meanHexNAcCoverage=0, percentUncovered=1, MS1_Score=0,
                    peptideLens=0, Obs_Mass=0, numStubs=-1):
    query_string = "(MS2_Score >= {MS2_Score} and  meanCoverage >= {meanCoverage} and\
     meanHexNAcCoverage >= {meanHexNAcCoverage} and \
     percentUncovered <= {percentUncovered} and  peptideLens >= {peptideLens}\
     and  Obs_Mass >= {Obs_Mass}\
     and  numStubs >= {numStubs})".format(**locals())
    return query_string


class PredictionResultsJSONEncoder(json.JSONEncoder):

    ''' '''

    def default(self, o):
        if isinstance(o, pd.DataFrame):
            return json.loads(o.to_json(orient="records"))
        else:
            return super(PredictionResultsJSONEncoder, self).default(o)


class IndirectIndexer(object):

    '''A stubby front-end to mimic the DataFrame indexer components like .ix[],
    .loc[], and such. Since these use __getitem__ and __setitem__ to do their
    work rather than direct function calls, must use an object'''

    def __init__(self, wrapped_indexer_getter, wrapped_indexer_setter):
        self.wrapped_indexer_getter = wrapped_indexer_getter
        self.wrapped_indexer_setter = wrapped_indexer_setter

    def __getitem__(self, inds):
        return self.wrapped_indexer_getter(inds)

    def __setitem__(self, inds, value):
        self.wrapped_indexer_setter(inds, value)


class PredictionResults(object):

    '''An extension-by-composition of pandas.DataFrame, carrying a dictionary of
    arbitrary metadata as well as a DataFrame of arbitrary values. Attempts to cover
    most basic use-cases, forwarding all attribute requests not satisfied by the wrapping
    composite to the DataFrame transparently. Operations that subset the DataFrame are
    captured and rewrapped to propagate metadata. Currently does not support merges or
    pandas's generic combination functions which do not pass through this class.'''

    __schema_version = 0x0001

    plugins = []

    def __init__(self, metadata, predictions):
        self.metadata = metadata
        self.predictions = predictions
        # Stub Indexers
        self.ix = IndirectIndexer(self.wrap_ix_get, self.wrap_ix_set)
        self.loc = IndirectIndexer(self.wrap_loc_get, self.wrap_loc_set)
        self.iloc = IndirectIndexer(self.wrap_iloc_get, self.wrap_iloc_set)

        version = self.metadata.get("schema_version", None)
        if version is None:
            self.schema_upgrade(version)
            self.metadata["schema_version"] = PredictionResults.__schema_version
            # May need to run some update code in the future
        if version != PredictionResults.__schema_version:
            self.schema_upgrade(version)
            # May need to run some update code in the future

    def schema_upgrade(self, version):
        pass

    def __getitem__(self, inds):
        val = self.predictions[inds]
        if len(self.index) == len(inds):
            val = PredictionResults(self.metadata, val)
        return val

    def __setitem__(self, inds, value):
        self.predictions[inds] = value

    def __getattr__(self, name):
        try:
            return object.__getattribute__(self, name)
        except:
            try:
                out = getattr(object.__getattribute__(self, "predictions"), name)
                if isinstance(out, MethodType):
                    out = partial(self._framewrapper_fn, name)
                return out
            except:
                raise AttributeError(
                    "Prediction Results has no attribute {attr}".format(attr=name))

    def __setattr__(self, name, value):
        if hasattr(self, "predictions"):
            if name in self:
                self[name] = value
            else:
                object.__setattr__(self, name, value)
        else:
            object.__setattr__(self, name, value)

    def __repr__(self):
        rep = "<PredictionResults>\n{frame}".format(frame=repr(self.predictions))
        return rep

    def __iter__(self):
        return iter(self.predictions)

    def __dir__(self):
        return ["metadata"] + [attr for attr in dir(self.predictions) if "__" not in attr]

    def __contains__(self, name):
        return name in self.predictions

    def wrap_ix_get(self, index):
        df = self.predictions.ix[index]
        return PredictionResults(self.metadata, df)

    def wrap_ix_set(self, index, value):
        self.predictions.ix[index] = value

    def wrap_loc_get(self, index):
        df = self.predictions.loc[index]
        return PredictionResults(self.metadata, df)

    def wrap_loc_set(self, index, value):
        self.predictions.loc[index] = value

    def wrap_iloc_get(self, index):
        df = self.predictions.iloc[index]
        return PredictionResults(self.metadata, df)

    def wrap_iloc_set(self, index, value):
        self.predictions.iloc[index] = value

    @property
    def fdr(self):
        return self.metadata.get("fdr", None)

    def fragments(self):
        return self.b_ions_with_HexNAc + self.y_ions_with_HexNAc + self.b_ion_coverage + self.y_ion_coverage

    def _framewrapper_fn(self, attr, *args, **kwargs):
        out = getattr(self.predictions, attr)(*args, **kwargs)
        if isinstance(out, pd.DataFrame):
            out = PredictionResults(self.metadata, out)
        return out

    def score_naive(self, peptide_weight=0.5, hexnac_weight=0.5,
                    uncovered_penalty_weight=1.0, stub_ion_weight=0.2,
                    backbone_fragment_weight=0.8):
        if (peptide_weight + hexnac_weight) != 1.0:
            raise ValueError("Peptide and HexNAc Weights must sum to 1.0")
        backbone_fragment_score = (self.meanCoverage * peptide_weight) +\
            (self.meanHexNAcCoverage * hexnac_weight) -\
            (uncovered_penalty_weight * self.percentUncovered)
        score = (backbone_fragment_score * backbone_fragment_weight) +\
            (stub_ion_weight * self.stubsObservedVsExpected)
        return score

    def optimize_fdr(self, n=1):
        if "fdr" not in self.metadata:
            raise KeyError("The FDR has not been calculated")
        fdr = self.metadata.get("fdr")
        # return fdr.ix[fdr.query("false_discovery_rate <= 0.05").num_real_matches.idxmax()]
        return fdr.query("false_discovery_rate <= 0.05").sort(
            ["num_real_matches"], ascending=False).iloc[:n]

    def kvquery(self, frame=None):
        if frame is None:
            frame = self.optimize_fdr()
        args = {k: v for k, v in frame.to_dict(orient="records")[0].items()
                if k not in {"false_discovery_rate", "num_real_matches", "num_decoy_matches"}}
        print(args)
        return self.query(query_threshold(**args))

    def serialize(self, target_buffer=None, close=True):
        if target_buffer is None:
            target_buffer = StringIO()
        if isinstance(target_buffer, basestring):
            target_buffer = open(target_buffer, "wb")
        encoder = PredictionResultsJSONEncoder()

        serialize_data = {
            "metadata": self.metadata,
            "predictions": self.predictions
        }

        for plugin in PredictionResults.plugins:
            plugin.serialize(self, serialize_data)

        for chunk in encoder.iterencode(serialize_data):
            target_buffer.write(chunk)
        if close:
            target_buffer.close()
        return

    @classmethod
    def deserialize(cls, data_buffer, recompute=False):
        if isinstance(data_buffer, basestring):
            data_buffer = open(data_buffer, "rb")
        loose_data = json.load(data_buffer)
        metadata = loose_data["metadata"]
        # if "fdr" in metadata:
        #     metadata["fdr"] = pd.DataFrame(metadata["fdr"])
        predictions = loose_data["predictions"]
        predictions = cls.ensure_fields(predictions, recompute=recompute)
        instance = cls(metadata, predictions)

        for plugin in PredictionResults.plugins:
            plugin.deserialize(instance, loose_data)

        return instance

    @classmethod
    def prepare(cls, metadata, predictions):
        prediction_frame = PredictionResults(metadata, predictions)
        prediction_frame.predictions = prediction_frame.ensure_fields(
            predictions)
        return prediction_frame

    @staticmethod
    def ensure_fields(raw_predictions, recompute=False):
        data = pd.DataFrame(raw_predictions)
        data.Peptide_mod = data.Peptide_mod.astype(str).replace("nan", "")
        # Handle Decoys
        data.ix[np.isnan(data.vol)].vol = -1
        try:
            data.ix[data.startAA.isnull()].startAA = -1
            data.ix[np.isnan(data.startAA)].startAA = -1
        except TypeError:
            data.startAA = -1
        try:
            data.ix[data.endAA.isnull()].endAA = -1
            data.ix[np.isnan(data.endAA)].endAA = -1
        except TypeError:
            data.endAA = -1
        if "bad_oxonium_ions" not in data or recompute:
            data["bad_oxonium_ions"] = data.Oxonium_ions.apply(lambda x: [])
        if "peptideLens" not in data or recompute:
            data = get_sequence_length(data)
        if "numOxIons" not in data or recompute:
            data["numOxIons"] = map(len, data["Oxonium_ions"])
        if "numStubs" not in data or recompute:
            data["numStubs"] = map(len, data["Stub_ions"])
        if "percent_b_ion_coverage" not in data or recompute:
            shim_percent_calcs(data)
        if "glycosylation_sites" not in data or recompute:
            data["glycosylation_sites"] = data.Glycopeptide_identifier.apply(find_occupied_glycosylation_sites)
        if "meanCoverage" not in data or recompute:
            ion_coverage_score(data)
        if "modificationSignature" not in data or recompute:
            modification_signature(data)
        if "glycanMass" not in data or recompute:
            infer_glycan_mass(data)
        if "stubsObservedVsExpected" not in data or recompute:
            data["stubsObservedVsExpected"] = data.apply(stubs_observed_expected_ratio, 1)
        return data


class PredictionResultsPlugin(object):

    @classmethod
    def serialize(self, instance, json_dict):
        pass

    @classmethod
    def deserialize(self, instance, json_dict):
        pass


class PredictionResultsPluginFDRPlugin(PredictionResultsPlugin):

    @classmethod
    def deserialize(self, instance, json_dict):
        if "fdr" in instance.metadata:
            instance.metadata["fdr"] = pd.DataFrame(instance.metadata["fdr"])
    # Does not need deserialize method, PredictionResultsJSONEncoder handles it
    # implicitly.
PredictionResults.plugins.append(PredictionResultsPluginFDRPlugin)


class PredictionResultsIntactMassNoFragmentsPlugin(PredictionResultsPlugin):

    @classmethod
    def deserialize(self, instance, json_dict):
        if "intact_mass_no_fragments" in json_dict:
            instance.intact_mass_no_fragments = json_dict["intact_mass_no_fragments"]

    @classmethod
    def serialize(self, instance, json_dict):
        try:
            json_dict["intact_mass_no_fragments"] = instance.intact_mass_no_fragments
        except:
            json_dict["intact_mass_no_fragments"] = []
PredictionResults.plugins.append(PredictionResultsIntactMassNoFragmentsPlugin)


class PredictionResultsSqlite(sqlitedict.SqliteDict):

    def __init__(self, file_name, frame):
        super(PredictionResultsSqlite, self).__init__(file_name)
        self.column_to_type = frame.apply(pd.lib.infer_dtype, 0).to_dict()
        for i, row in frame.iterrows():
            self[i] = row.to_dict()

    def __getitem__(self, key):
        val = super(PredictionResultsSqlite, self).__getitem__(key)
        return pd.Series(val)

    def __iter__(self):
        for i, row in self.items():
            yield pd.Series(row)


def convert_csv_to_nested(path):
    frame = prepare_model_file(path)
    metadata = {
        # Shim in the ms1_results_file path
        "ms1_results_file": path,
        "_converted": True
    }
    pr = PredictionResults(metadata, frame)
    return pr
