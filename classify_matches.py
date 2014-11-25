from os.path import splitext
import numpy as np
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression

from statistics import *


class ModelTask(object):

    '''Base Class for representing tasks ran about a model.'''
    method_table = methods

    @classmethod
    def dump_options(cls):
        return cls.method_table.keys()

    def __init__(self, method='full_random_forest', method_init_args=None, method_fit_args=None):
        if method_init_args is None:
            method_init_args = dict()
        if method_fit_args is None:
            method_fit_args = dict()
        self.method = method
        self.method_init_args = method_init_args
        self.method_fit_args = method_fit_args


class PrepareModelTask(ModelTask):

    def __init__(self, model_file_path, output_path=None, method="full_random_forest",
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
        self.model_frame['noise_filter'] = 0.0

    def save_model(self):
        save_model_file(self.model_frame, self.output_path)

    def run(self):
        self.build_model()
        self.save_model()
        return self.output_path


class ClassifyTargetWithModelTask(ModelTask):

    def __init__(
        self, model_file_path, target_file_path, output_path=None, method="full_random_forest",
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
        self.target_frame["noise_filter"] = scores
        self.target_frame["MS2_Score"] = scores * (((self.target_frame.meanCoverage +
                                                     self.target_frame.meanHexNAcCoverage) / 2.0)
                                                   - self.target_frame.percentUncovered)
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
        self, model_file_path, method="full_random_forest",
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
            ax.text(
                0.5, 0.5, "This classifier does not support\nfeature importance")
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
                               [0, 1], model_conf_matrix[1, 0], model_conf_matrix[1, 1]]
        null_conf_measures = [null_conf_matrix[0, 0], null_conf_matrix
                              [0, 1], null_conf_matrix[1, 0], null_conf_matrix[1, 1]]

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


class CompareModelsDiagnosticTask(ModelDiagnosticsTask):

    def __init__(self, model_file_path, labeled_file_path, method="full_random_forest",
                 output_path=None, method_init_args=None, method_fit_args=None):
        if output_path is None:
            output_path = splitext(
                splitext(model_file_path)[0])[0] + ".errors_%s.png" % method
        super(CompareModelsDiagnosticTask, self).__init__(model_file_path, method,
                                                          output_path, method_init_args, method_fit_args)
        self.labeled_file_path = labeled_file_path

    def plot_error_rate(self, ax):
        width = 0.3
        labeled_frame = prepare_model_file(self.labeled_file_path)
        model_conf_matrix = generate_confusion_matrix(
            self.classifier, labeled_frame, self.model_formula)
        # True Positive, False Negative, False Positive, True Negative
        model_conf_measures = [model_conf_matrix[0, 0], model_conf_matrix
                               [0, 1], model_conf_matrix[1, 0], model_conf_matrix[1, 1]]

        #print((self.model_frame['call'] == True).sum(), self.model_frame['call'].count())

        model_rects = ax.bar(np.arange(len(model_conf_measures)),
                             model_conf_measures,
                             width=width, align="center", color='green')

        ax.set_title("Classifier Performance Rates")
        ax.set_xticks(np.arange(len(model_conf_measures)))
        ax.set_xticklabels(["TPR", "FPR", "FNR", "TNR"])
        ax.set_ylabel("Number of observations")
        ax.autoscale()
