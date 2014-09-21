
#import io
import re
from os.path import splitext
import json

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.metrics import roc_curve, auc, confusion_matrix
from sklearn.cross_validation import ShuffleSplit

# Defines the list of fields that must be JSON encoded and decoded
# for storage in CSV files. Because regular scientists cannot possibly
# understand the fractal madness of pure hierarchical structured data
# files. YAML is out of the question.
json_serialized_columns = ["Oxonium_ions", "Stub_ions",
                           "b_ion_coverage", "b_ions_with_HexNAc",
                           "y_ion_coverage", "y_ions_with_HexNAc",
                           "peptideCoverageMap", "hexNAcCoverageMap"]

# A map for translating legacy labels for 'call' to and from
# Python booleans
call_map = {
    "Yes": True,
    "No": False,
    True: "Yes",
    False: "No"
}


# Decode :json_serialized_columns columns into Python objects for interaction.
# These columns are not "fast", and should only be accessed for
# big tasks. This function destructively modifies the input object.
def deserialize_compound_fields(data_struct):
    for field in json_serialized_columns:
        try:
            data_struct[field] = map(json.loads, data_struct[field])
        except TypeError:
            pass
        except KeyError:
            pass

    return data_struct


# Encode :json_serialized_columns columns into JSON strings for writing
# to CSV. This is a destructive operation and should be called on a copy
# of the working object, or else the object will need to go through
# :deserialize_compound_fields again. This function destructively modifies
# the input object.
def serialize_compound_fields(data_struct):
    for field in json_serialized_columns:
        data_struct[field] = map(json.dumps, data_struct[field])
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


# Compute the peptide sequence lengths, which is harder than taking the len()
# of the sequence string. Must exclude all modifications in "()" and the glycan
# composition at the end with "[]".
def get_sequence_length(data):
    data["peptideLens"] = data.Peptide.str.split(
        r'\(.*?\)').str.join(' ').str.len()
    return data


# Given the path to a CSV file for a model or target to be classified, parse the
# CSV into a pandas.DataFrame and make sure it has all of the necessary columns and
# transformations to either build a model with it or
def prepare_model_file(path):
    data = pd.read_csv(path)
    deserialize_compound_fields(data)
    data.Peptide_mod = data.Peptide_mod.astype(str).str.replace("nan", "")
    data['abs_ppm_error'] = data.ppm_error.abs()
    if "peptideLens" not in data:
        data = get_sequence_length(data)
    if "numOxIons" not in data:
        data["numOxIons"] = map(len, data["Oxonium_ions"])
    if "numStubs" not in data:
        data["numStubs"] = map(len, data["Stub_ions"])
    if "percent_b_ion_coverage" not in data:
        shim_percent_calcs(data)
    get_ion_coverage_score(data)

    return data


# Serialize the model object into a CSV file. Creates a copy of the data which has its
# complex fields JSON encoded.
def save_model_file(data_struct, path):
    write_copy = data_struct.copy()
    serialize_compound_fields(write_copy)
    write_copy.to_csv(path, index=False, encoding="utf-8")


def get_ion_coverage_score(data_struct):
    coverage_data = data_struct.apply(compute_ion_coverage_map, 1)
    data_struct['meanCoverage'] = coverage_data["meanCoverage"]
    data_struct['percentUncovered'] = coverage_data["percentUncovered"]
    data_struct['meanHexNAcCoverage'] = coverage_data['meanHexNAcCoverage']
    data_struct['peptideCoverageMap'] = coverage_data['peptideCoverageMap']
    data_struct['hexNAcCoverageMap'] = coverage_data['hexNAcCoverageMap']
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
        y_ion_coverage += ion

    total_coverage += b_ion_coverage + b_ions_with_HexNAc
    total_coverage += y_ion_coverage[
        ::-1] + y_ions_with_HexNAc[::-1]  # Reverse

    percent_uncovered = (sum(total_coverage == 0) / float(peptide_length))
    mean_coverage = total_coverage.mean() / float(peptide_length)

    hexnac_coverage = (b_ions_with_HexNAc + y_ions_with_HexNAc[::-1])
    mean_hexnac_coverage = hexnac_coverage.mean() / float(peptide_length)

    return pd.Series({"meanCoverage": mean_coverage,
                     "percentUncovered": percent_uncovered,
                     "meanHexNAcCoverage": mean_hexnac_coverage,
                     "peptideCoverageMap": list(total_coverage),
                     "hexNAcCoverageMap": list(hexnac_coverage)})


# Based upon the assumption that the same MS1 Score + Observed Mass implies
# the same precursor ion, the data is grouped by these two scores, and sort within groups.
# If the input DataFrame does not have an MS2 Score column, it will be removed
# from the within-group sorting criteria.
def determine_ambiguity(data_struct):
    data_struct['ambiguity'] = False
    groupings = data_struct.groupby(["MS1_Score", "Obs_Mass"])
    # Bit list for pandas.DataFrame sort function
    sort_direction = [1, 1, 1, 0, 1, 0, 1]
    sort_fields = ["MS2_Score", "numStubs", "meanCoverage",
                   "percentUncovered", "peptideLens",
                   "abs_ppm_error", "meanHexNAcCoverage"]
    if "MS2_Score" not in data_struct:
        sort_direction.pop(0)
        sort_fields.pop(0)
    regrouping = []
    for name, group in groupings:
        if (len(group.index) == 1):
            regrouping.append(group)
            continue
        group.sort(columns=sort_fields, ascending=sort_direction,
                   axis=0, inplace=True)
        group['ambiguity'] = True
        regrouping.append(group)

    regroup = pd.concat(regrouping)
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

model_definitions = {
    "full": [
        "meanCoverage",
        "percentUncovered",
        "abs_ppm_error",
        "meanHexNAcCoverage",
        "peptideLens"
    ],
    "backbone_hexnac": [
        "meanHexNAcCoverage",
        "percentUncovered"
    ],
}


def fit_random_forest(
    data_struct, formula=None, max_features='auto', oob_score=True,
    n_estimators=100, criterion="entropy", n_jobs=2, model_init_args=None,
        model_fit_args=None):
    if formula is None:
        formula = model_definitions['full']
    if model_init_args is None:
        model_init_args = dict()
    if model_fit_args is None:
        model_fit_args = dict()
    model_init_args['max_features'] = max_features
    model_init_args['oob_score'] = oob_score
    model_init_args['n_jobs'] = 2
    model_init_args['n_estimators'] = n_estimators
    rfc = RandomForestClassifier(**model_init_args)
    rfc.fit(data_struct[formula], data_struct['call'], **model_fit_args)

    return rfc


def fit_gradient_boosted_forest(data_struct, formula=None, max_features="auto",
                                model_init_args=None, model_fit_args=None):
    if formula is None:
        formula = model_definitions['full']
    if max_features is None:
        max_features = len(formula)
    if model_init_args is None:
        model_init_args = dict()
    if model_fit_args is None:
        model_fit_args = dict()
    model_init_args['max_features'] = max_features
    gbc = GradientBoostingClassifier(**model_fit_args)
    gbc.fit(data_struct[formula], data_struct['call'], **model_fit_args)

    return gbc


def fit_radial_basis_svm(data_struct, formula=None, model_init_args=None, model_fit_args=None):
    if formula is None:
        formula = model_definitions['full']
    if model_init_args is None:
        model_init_args = dict()
    if model_fit_args is None:
        model_fit_args = dict()
    model_init_args['probability'] = True
    svc = SVC(**model_init_args)
    svc.fit(data_struct[formula], data_struct['call'], **model_fit_args)
    return svc


def classify_with_model(classifier, data_struct, formula=None):
    if formula is None:
        formula = model_definitions['full']
    class_probabilities = classifier.predict_proba(data_struct[
                                                   formula])
    return class_probabilities.T[1]


def generate_confusion_matrix(classifier, data_struct, formula=None, labels=None):
    if formula is None:
        formula = model_definitions['full']
    if labels is None and "call" not in data_struct and "MS2_Score" not in data_struct:
        raise Exception(
            "Data not labeled or scored. Cannot generate Confusion Matrix")
    if labels is None and "call" not in data_struct and "MS2_Score" in data_struct:
        labels = data_struct['MS2_Score'] > 0.5
    elif labels is None:
        labels = data_struct['call']
    pred = classifier.predict(data_struct[formula])
    confusion = confusion_matrix(labels, pred, [True, False])
    return confusion


def generate_null_model(classifier, data_struct, formula=None, model_fit_args=None):
    if formula is None:
        formula = model_definitions["full"]
    if model_fit_args is None:
        model_fit_args = dict()

    null_classifier = clone(classifier)
    random_labels = np.random.uniform(0, 1, len(data_struct)) > .5
    null_classifier.fit(data_struct[formula], random_labels, **model_fit_args)
    return null_classifier, random_labels


class ModelTask(object):

    '''Base Class for representing tasks ran about a model.'''
    method_table = {
        "full": ["full", fit_random_forest],
        "gradient_boosted": ["full", fit_gradient_boosted_forest],
        "backbone_hexnac": ["backbone_hexnac", fit_random_forest],
        "radial_basis_svm": ["full", fit_radial_basis_svm]
    }

    @classmethod
    def dump_options(cls):
        return cls.method_table.keys()

    def __init__(self, method='full', method_init_args=None, method_fit_args=None):
        if method_init_args is None:
            method_init_args = dict()
        if method_fit_args is None:
            method_fit_args = dict()
        self.method = method
        self.method_init_args = method_init_args
        self.method_fit_args = method_fit_args


class PrepareModelTask(ModelTask):

    def __init__(self, model_file_path, output_path=None, method="full",
                 method_init_args=None, method_fit_args=None):
        if method_init_args is None:
            method_init_args = dict()
        if method_fit_args is None:
            method_fit_args = dict()
        if output_path is None:
            output_path = splitext(
                splitext(model_file_path)[0])[0] + ".model.csv"
        self.output_path = output_path
        super(PrepareModelTask, self).__init__(method)
        self.model_path = model_file_path
        self.model_frame = prepare_model_file(model_file_path)

    def build_model(self):
        calls = call_by_coverage(self.model_frame)
        self.model_frame['call'] = calls
        self.model_frame = determine_ambiguity(self.model_frame)
        self.model_frame['MS2_Score'] = 0.0

    def save_model(self):
        save_model_file(self.model_frame, self.output_path)

    def run(self):
        self.build_model()
        self.save_model()
        return self.output_path


class ClassifyTargetWithModelTask(ModelTask):

    def __init__(
        self, model_file_path, target_file_path, output_path=None, method="full",
            method_init_args=None, method_fit_args=None):
        if method_init_args is None:
            method_init_args = dict()
        if method_fit_args is None:
            method_fit_args = dict()
        if output_path is None:
            output_path = splitext(
                splitext(target_file_path)[0])[0] + ".scored.csv"
        self.output_path = output_path
        super(ClassifyTargetWithModelTask, self).__init__(
            method, method_init_args, method_fit_args)
        self.model_path = model_file_path
        self.model_frame = prepare_model_file(model_file_path)
        self.target_path = target_file_path
        self.target_frame = prepare_model_file(target_file_path)
        self.model_formula, self.classifier_fn = ClassifyTargetWithModelTask.method_table[
            self.method]
        self.model_formula = model_definitions[self.model_formula]
        self.classifier = None

    def build_model(self):
        self.classifier = self.classifier_fn(
            self.model_frame, self.model_formula,
            model_init_args=self.method_init_args,
            model_fit_args=self.method_fit_args)

    def classify_with_model(self):
        scores = classify_with_model(
            self.classifier, self.target_frame, self.model_formula)
        self.target_frame["MS2_Score"] = scores
        self.target_frame = determine_ambiguity(self.target_frame)
        self.target_frame['call'] = scores > 0.5

    def save_target(self):
        save_model_file(self.target_frame, self.output_path)

    def run(self):
        self.build_model()
        self.classify_with_model()
        self.save_target()
        return self.output_path


class ModelDiagnosticsTask(ModelTask):

    def __init__(
        self, model_file_path, method="full",
            output_path=None,
            method_init_args=None, method_fit_args=None):
        if method_init_args is None:
            method_init_args = dict()
        if method_fit_args is None:
            method_fit_args = dict()
        if output_path is None:
            output_path = splitext(
                splitext(model_file_path)[0])[0] + ".feature_importance_%s.png" % method
        self.output_path = output_path
        super(ModelDiagnosticsTask, self).__init__(method)
        self.model_path = model_file_path
        self.model_frame = prepare_model_file(model_file_path)
        self.model_formula, self.classifier_fn = ModelDiagnosticsTask.method_table[
            method]
        self.model_formula = model_definitions[self.model_formula]
        self.classifier = None

    def build_model(self):
        self.classifier = self.classifier_fn(
            self.model_frame, self.model_formula,
            model_init_args=self.method_init_args,
            model_fit_args=self.method_fit_args)

    def prepare_figure(self):
        plt.figure(figsize=(6, 12))
        fig, (ax) = plt.subplots(1, 2)
        self.axes = ax

    def plot_feature_importance(self, ax):
        if not hasattr(self.classifier, "feature_importances_"):
            #raise Exception("Has no feature importances")
            ax.text(0.5, 0.5, "This classifier does not support\nfeature importance")
            ax.axis([0, 1, 0, 1])
        else:
            ax.bar(np.arange(len(self.classifier.feature_importances_)),
                   self.classifier.feature_importances_, align='center')
            ax.set_ylabel("Feature Importance")
            ax.set_xticks(np.arange(len(self.classifier.feature_importances_)))
            ax.set_xticklabels(self.model_formula, rotation=45, ha='right')
            ax.set_title("Classifier Feature Importance")

    def plot_error_rate(self, ax):
        width = 0.3

        model_conf_matrix = generate_confusion_matrix(
            self.classifier, self.model_frame, self.model_formula)
        null_classifier, random_labels = generate_null_model(
            self.classifier, self.model_frame, self.model_formula, self.method_fit_args)
        assert not all(np.logical_and(random_labels, self.model_frame['call']))
        null_conf_matrix = generate_confusion_matrix(
            self.classifier, self.model_frame, self.model_formula, random_labels)

        # True Positive, False Negative, False Positive, True Negative
        model_conf_measures = [model_conf_matrix[0, 0], model_conf_matrix
                              [0, 1], model_conf_matrix[1, 0], model_conf_matrix[0, 1], model_conf_matrix[1, 1]]
        null_conf_measures = [null_conf_matrix[0, 0], null_conf_matrix
                             [0, 1], null_conf_matrix[1, 0], null_conf_matrix[0, 1], null_conf_matrix[1, 1]]

        #print((self.model_frame['call'] == True).sum(), self.model_frame['call'].count())

        model_rects = ax.bar(np.arange(len(model_conf_measures)),
                             model_conf_measures,
                             width=width, align="center", color='green')
        null_rects = ax.bar(np.arange(len(null_conf_measures)) + width,
                            null_conf_measures, width=width, align="center", color='red')

        ax.set_title("Classifier Performance Rates")
        ax.set_xticks(np.arange(len(model_conf_measures)) + width)
        ax.set_xticklabels(["TPR", "FPR", "FNR", "TNR"])
        ax.set_ylabel("Number of observations")
        ax.legend((model_rects[0], null_rects[0]), ("Model", "Null"))
        ax.autoscale()

    def save_figure(self):
        plt.tight_layout()
        plt.savefig(self.output_path)
        pass

    def run(self):
        self.build_model()
        self.prepare_figure()
        self.plot_feature_importance(self.axes[0])
        self.plot_error_rate(self.axes[1])
        self.save_figure()
        return self.output_path
