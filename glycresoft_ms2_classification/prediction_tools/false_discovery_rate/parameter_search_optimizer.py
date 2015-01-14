import itertools
import collections
import logging
import operator

logger = logging.getLogger(__name__)
import numpy as np
import pandas as pd

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
    "peptideLens": int_search(inc(3)),
    "percentUncovered": float_search(inc(0.05)),
    "meanCoverage": float_search(inc(0.05)),
    "meanHexNAcCoverage": float_search(inc(0.05)),
}


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


def calculate_fdr_task(paramquery, predictions, decoys, decoy_ratio, min_discoveries=10):
    param, query = paramquery
    print(param)
    n_predictions = len(predictions.query(query).index)
    n_decoys = len(decoys.query(query).index)
    if n_decoys < 1:
        fdr = 0
    if n_predictions <= min_discoveries:
        return (0, 1)
    else:
        fdr = (n_decoys/float(n_predictions + n_decoys)) * (1 + (decoy_ratio))
    results_obj = {
        "num_real_matches": n_predictions,
        "num_decoy_matches": n_decoys,
        "decoy_to_real_ratio": decoy_ratio,
        "false_discovery_rate": fdr,
        "threshold": query,
    }
    try:
        results_obj.update(param)
    except:
        pass
    print("Result:", str(results_obj))
    return (str(results_obj))


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


class GridSearchOptimizer(object):
    def __init__(self, predictions, decoys, decoy_ratio=20.0, filter_terms=None, conditions=None):
        self.trials = {}
        if filter_terms is None:
            filter_terms = list(filter_terms_type_map.keys())
        if conditions is None:
            conditions = {}
        self.predictions = predictions
        self.decoys = decoys
        self.terms = [self.make_search_param(t, filter_terms_type_map.get(t, t)) for t in filter_terms]
        self.conditions = conditions
        self.Parameters = collections.namedtuple("Parameters", [t['name'] for t in self.terms])
        self.decoy_ratio = float(decoy_ratio)

    def make_search_param(self, name, strategy):
        if isinstance(name, dict):
            name = dict.get("name")
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
        total = reduce(map(len, ranges[1:]), operator.pow, len(ranges[0])) if len(ranges) > 1 else len(ranges[0])
        i = 0
        logger.info("Parameter Space contains %d combinations", total)
        for combination in itertools.product(*ranges):
            i += 1
            yield self.Parameters(*combination)
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
        self.trials[param] = results_obj
        return fdr

    def run_all(self, min_discoveries=10):
        print("start runall")
        for param in self.param_space():
            self.calculate_fdr(param, min_discoveries)
        print("end runall")

    def select_best_criterion(self, n=5):
        trials = pd.DataFrame(self.trials.values())
        print(trials.columns)
        self._trials = trials.sort(columns=["false_discovery_rate", "num_real_matches"], ascending=[1, 0])
        return self._trials.iloc[:n]

    def optimize(self, n=5, min_discoveries=10):
        self.run_all(min_discoveries)
        return self.select_best_criterion(n)

    def gradient_descent(self, convergence_delta=0.001, burn_in=1000):
        best_score = 1.0
        best_scores = [100, 200, 300, 400] * 3
        param = self.Parameters(*[t['start'] for t in self.terms])
        print(param)
        def convergence():
            #print(best_scores)
            return all((b - best_score) <= convergence_delta for b in best_scores)

        cnt = 0
        while(not (convergence() and cnt > burn_in)):
            cnt += 1
            new_param, new_score = gradient_descent(list(step_param(param, self.terms)), self.calculate_fdr)
            if new_score < best_scores:
                next_ = best_scores[0]
                for i in range(len(best_scores) - 1):
                    asn = best_scores[i + 1]
                    best_scores[i + 1] = next_
                    next_ = asn
                best_scores[0] = best_score
                best_score = new_score
                param = new_param
                print("New Parameters: {0}, New Score: {1}".format(new_param, best_score))


def grid_search(predictions, decoys, filter_terms=None):
    pass


def gradient_descent(parameters, scoring_function, selector=min):
    scores = [(param, scoring_function(param)) for param in parameters if scoring_function(param) is not None]
    best_param, score = selector(scores, key=lambda x: x[1]) if len(scores) > 0 else None, float('inf')
    return best_param, score
