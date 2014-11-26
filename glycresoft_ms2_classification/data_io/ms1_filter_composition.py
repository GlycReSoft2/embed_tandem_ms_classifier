import argparse
import os

import pandas as pd


def filter_composition_hypothesis_by_ms1_results(composition_hypothesis, ms1_results, score_threshold):
    composition_groups = composition_hypothesis.groupby(["Compositions", "Peptide Sequence", "Peptide Modification"])
    ms1_groups = ms1_results.groupby(["Compound Key", "PeptideSequence", "PeptideModification"])
    accepted_indices = []
    keys = set(ms1_groups.groups.keys()) & set(composition_groups.groups.keys())
    for key in keys:
        results = ms1_groups.get_group(key)
        if any(results.Score >= score_threshold):
            try:
                accepted_indices.extend(composition_groups.groups[key])
            except Exception, e:
                print(e)
                print(results)

    return composition_hypothesis.ix[accepted_indices]


def main(args):
    composition_hypothesis = pd.read_csv(args.composition_hypothesis_path)
    composition_hypothesis["Peptide Sequence"] = composition_hypothesis["Peptide Sequence"].astype(str).replace("nan", "")
    composition_hypothesis["Compositions"] = composition_hypothesis["Compositions"].astype(str).replace("nan", "")
    composition_hypothesis["Peptide Modification"] = composition_hypothesis["Peptide Modification"].astype(str).replace("nan", "")
    ms1_results = pd.read_csv(args.ms1_results_path)
    ms1_results.PeptideModification = ms1_results.PeptideModification.astype(str).replace("nan", "")
    ms1_results.PeptideSequence = ms1_results.PeptideSequence.astype(str).replace("nan", "")
    ms1_results["Compound Key"] = ms1_results["Compound Key"].astype(str).replace("nan", "")
    threshold = args.threshold
    filtered = filter_composition_hypothesis_by_ms1_results(composition_hypothesis, ms1_results, threshold)
    if args.out is None:
        args.out = os.path.splitext(args.composition_hypothesis_path)[0] + ".filtered.csv"
    filtered.to_csv(args.out, index=False)


app = argparse.ArgumentParser()
app.add_argument("-c", "--composition_hypothesis_path")
app.add_argument("-m", "--ms1_results_path")
app.add_argument("-s", "--threshold", default=0.5, type=float)
app.add_argument("-o", "--out", required=False, default=None)

if __name__ == '__main__':
    main(app.parse_args())
