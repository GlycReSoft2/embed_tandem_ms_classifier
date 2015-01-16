import itertools
from collections import OrderedDict

import numpy as np

from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix

model_definitions = {
    "full": [
        "meanCoverage",
        "percentUncovered",
        "meanHexNAcCoverage",
        "peptideLens",
        "stubsObservedVsExpected"
    ],
    "backbone_hexnac": [
        "meanHexNAcCoverage",
        "percentUncovered",
    ],
}

fitting_functions = []


def fit_random_forest(
    data_struct, formula=None, max_features='auto', oob_score=True,
    n_estimators=100, criterion="entropy", n_jobs=1, model_init_args=None,
        model_fit_args=None):
    if formula is None:
        formula = model_definitions['full']
    if model_init_args is None:
        model_init_args = dict()
    if model_fit_args is None:
        model_fit_args = dict()
    model_init_args['max_features'] = max_features
    model_init_args['oob_score'] = oob_score
    model_init_args['n_jobs'] = n_jobs
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


def fit_svc(data_struct, formula=None, model_init_args=None, model_fit_args=None):
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


def fit_logistic_regression(data_struct, formula=None, model_init_args=None, model_fit_args=None):
    if formula is None:
        formula = model_definitions['full']
    if model_init_args is None:
        model_init_args = dict()
    if model_fit_args is None:
        model_fit_args = dict()
    model_init_args["penalty"] = model_init_args.get("penalty", "l2")
    log_reg = LogisticRegression(**model_init_args)
    log_reg.fit(data_struct[formula], data_struct["call"], **model_fit_args)
    return log_reg

fitting_functions.extend(
    [fit_random_forest, fit_gradient_boosted_forest, fit_svc, fit_logistic_regression])


methods = OrderedDict()
for model_name, fitter in itertools.product(model_definitions.keys(), fitting_functions):
    methods["%s_%s" % (model_name, fitter.func_name.replace("fit_", ""))] = (model_name, fitter)
methods['naive'] = ("naive", "naive")


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
    random_labels = np.random.uniform(0, 1, len(data_struct.index)) > .5
    null_classifier.fit(data_struct[formula], random_labels, **model_fit_args)
    return null_classifier, random_labels
