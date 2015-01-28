import itertools
import logging
import operator

logger = logging.getLogger(__name__)
import numpy as np
import pandas as pd
import xray
# Search Strategies

int_search = lambda s: {"type": int, "start": "min", "step": s, "stop": "max"}
float_search = lambda s: {"type": float, "start": "min", "stop": "max", "step": s}

gte = " >= "
lt = " < "

inc = lambda s: lambda x: x + s

filter_terms_type_map = {
    "MS2_Score": float_search("stream"),
    "MS1_Score": float_search("stream"),
    "ppm_error": float_search(inc(0.5)),
    "numStubs": int_search(inc(1)),
    "Obs_Mass": float_search(inc(100)),
    "peptideLens": int_search("stream"),
    "percentUncovered": float_search(inc(0.05)),
    "meanCoverage": float_search(inc(0.05)),
    "meanHexNAcCoverage": float_search(inc(0.05)),
}


def fdr1(num_real, num_decoys, *args, **kwargs):
    return num_decoys / float(num_real)

def fdr2(num_real, num_decoys, *args, **kwargs):
    return (num_decoys * 2)/float(num_real + num_decoys)

def fdr3(num_real, num_decoys, ratio, *args, **kwargs):
    return (num_decoys * (1 + (1/float(ratio)))) / float(num_real + num_decoys)


def make_search_strategy(name, start="min", step_fn=inc(0.05), stop="max"):
    strategy = {
        "name": name,
        "start": start,
        "step": step_fn,
        "stop": stop
    }
    return strategy


def query_threshold(MS2_Score=0, meanCoverage=0, meanHexNAcCoverage=0, percentUncovered=1, MS1_Score=0,
                    peptideLens=0, Obs_Mass=0, numStubs=-1):
    query_string = "(MS2_Score >= {MS2_Score} and  meanCoverage >= {meanCoverage} and\
     meanHexNAcCoverage >= {meanHexNAcCoverage} and \
     percentUncovered <= {percentUncovered} and  peptideLens >= {peptideLens}\
     and  Obs_Mass >= {Obs_Mass}\
     and  numStubs >= {numStubs})".format(**locals())
    return query_string


def step_param(param, terms):
    for term in terms:
        current = getattr(param, term["name"])
        current = term['step'](current)
        if current > term['stop'] or current < term['start']:
            continue
        yield param._replace(**{term["name"]: current})


def stream_steps(*data):
    dataset = (sorted(list(reduce(operator.or_, map(set, data), set()))))
    _iter = iter(dataset)

    def stream(*args, **kwargs):
        return _iter.next()

    return stream, dataset


class CountExclusion(object):
    def __init__(self, predictions, decoys, decoy_ratio=20, filter_terms=None):
        if filter_terms is None:
            filter_terms = list(filter_terms_type_map.values())
        self.predictions = predictions
        self.decoys = decoys
        self.terms = [self.make_search_param(t) for t in filter_terms]
        self.decoy_ratio = float(decoy_ratio)

    def traverse(self):
        ranges = [p['range'] for p in self.terms]
        names = [p['name'] for p in self.terms]
        dims = map(len, ranges)
        prediction_array = xray.DataArray(np.zeros(dims), dims=names, coords=ranges)
        decoy_array = xray.DataArray(np.zeros(dims), dims=names, coords=ranges)
        if(len(names) > 1):
            apply_fn = map
        else:
            apply_fn = lambda x, y: (x(y),)

        def slice_back(i):
            return slice(None, i)
        for feature_levels, owned_indices in self.predictions.groupby(names).groups.items():
            prediction_array.loc[tuple(apply_fn(slice_back, feature_levels))] += len(owned_indices)

        for feature_levels, owned_indices in self.decoys.groupby(names).groups.items():
            decoy_array.loc[tuple(apply_fn(slice_back, feature_levels))] += len(owned_indices)

        self.trials = np.divide(decoy_array * (1 + (1/self.decoy_ratio)), decoy_array + prediction_array).to_dataframe()
        self.trials.reset_index(drop=False, inplace=True)
        self.trials.rename(columns={None: "false_discovery_rate"}, inplace=True)
        self.trials["num_real_matches"] = np.ravel(prediction_array)
        self.trials["num_decoy_matches"] = np.ravel(decoy_array)
        self.trials.dropna(0, 'any', inplace=True)
        self.trials.sort(["false_discovery_rate", "num_real_matches"], ascending=[True, False], inplace=True)

    def make_search_param(self, name):
        param = {}
        param['step'], param['range'] = stream_steps(self.predictions[name], self.decoys[name])
        param['name'] = name
        return param

    def select_best_criterion(self, n=5, threshold=0.05):
        return self.trials.query(
            "false_discovery_rate <= {0}".format(threshold)).sort(
            ["num_real_matches"], ascending=False).iloc[:n]

    def compress(self, threshold=0.05):
        return pd.DataFrame(
            cases.sort("num_real_matches", ascending=False).iloc[0]
            for fdr, cases in self.trials.query(
                "false_discovery_rate <= {0}".format(threshold)).groupby("false_discovery_rate"))

    def optimize(self, n=5, threshold=0.05, min_discoveries=10):
        self.traverse()
        return self.select_best_criterion(n, threshold)


class GridSearchOptimizer(object):
    def __init__(self, predictions, decoys, decoy_ratio=20.0, filter_terms=None, conditions=None):
        self.trials = {}
        if filter_terms is None:
            filter_terms = list(filter_terms_type_map.values())
        if conditions is None:
            conditions = {}
        self.predictions = predictions
        self.decoys = decoys
        self.terms = [self.make_search_param(t) for t in filter_terms]
        self.conditions = conditions
        self.decoy_ratio = float(decoy_ratio)

    def make_search_param(self, strategy):
        if isinstance(strategy, basestring):
            strategy = filter_terms_type_map[strategy]
        name = strategy.get("name")
        param = dict(strategy)
        start = strategy["start"]
        if start == "max":
            start = max(self.predictions[name].max(), self.decoys[name].max())
        elif start == "min":
            start = min(self.predictions[name].min(), self.decoys[name].min())
        param["start"] = start
        stop = strategy["stop"]
        if stop == "max":
            stop = max(self.predictions[name].max(), self.decoys[name].max())
        elif stop == "min":
            stop = min(self.predictions[name].min(), self.decoys[name].min())
        param["stop"] = stop
        if param['step'] == "stream":
            param['step'], param['range'] = stream_steps(self.predictions[name], self.decoys[name])
        else:
            param['range'] = np.arange(param["start"], param["step"](param["stop"]), param["step"](0))
        param['name'] = name
        return param

    def param_space(self):
        ranges = [p['range'] for p in self.terms]
        names = [p['name'] for p in self.terms]
        total = reduce(operator.pow, map(len, ranges[1:]), len(ranges[0])) if len(ranges) > 1 else len(ranges[0])
        i = 0
        logger.info("Parameter Space contains %d combinations", total)
        for combination in itertools.product(*ranges):
            i += 1
            yield dict(zip(names, combination))
            if i % (total / 10) == 0:
                logger.info("Searched %f%% of parameter space", (i / float(total)) * 100)

    def make_query(self, params):
        param_dict = params._asdict() if not isinstance(params, dict) else params
        param_dict.update(self.conditions)
        return query_threshold(**param_dict)

    def calculate_fdr(self, param, min_discoveries=10):
        if param in self.trials:
            return self.trials[param]['false_discovery_rate']
        query = self.make_query(param)
        n_predictions = len(self.predictions.query(query).index)
        n_decoys = len(self.decoys.query(query).index)
        if n_decoys < 1:
            fdr = 0
        if n_predictions <= min_discoveries:
            return None
        else:
            fdr = (n_decoys/float(n_predictions + n_decoys)) * (1 + (1/self.decoy_ratio))
        results_obj = {
            "num_real_matches": n_predictions,
            "num_decoy_matches": n_decoys,
            "decoy_to_real_ratio": self.decoy_ratio,
            "false_discovery_rate": fdr,
            "threshold": query,
        }
        results_obj.update(param._asdict())
        self.trials[frozenset(param.items())] = results_obj
        return fdr

    def run_all(self, min_discoveries=10):
        print("start runall")
        for param in self.param_space():
            self.calculate_fdr(param, min_discoveries)
        print("end runall")

    def select_best_criterion(self, n=5, threshold=0.05):
        trials = pd.DataFrame(self.trials.values())
        self._trials = trials.sort(columns=["false_discovery_rate", "num_real_matches"], ascending=[1, 0])
        return self._trials.query("false_discovery_rate <= {0}".format(threshold)).sort(["num_real_matches"])

    def optimize(self, n=5, threshold=0.05, min_discoveries=10):
        self.run_all(min_discoveries)
        return self.select_best_criterion(n, threshold)
