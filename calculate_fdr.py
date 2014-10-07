import csv
import itertools
import random
import functools

from collections import defaultdict

import pandas as pd

from random_glycopeptide import random_glycopeptide_to_sequence_space
import classify_matches
import match_ions
import postprocess
from structure.sequence import Sequence


def make_predicate(**kwargs):
    return functools.partial(predicate_base, **kwargs)


def predicate_base(x, MS2_Score=0, meanCoverage=0, meanHexNAcCoverage=0, percentUncovered=1, MS1_Score=0, vol=-1,
                   peptideLens=0, Obs_Mass=0):
    return (x.MS2_Score >= MS2_Score) & (x.meanCoverage >= meanCoverage) & (x.percentUncovered <= percentUncovered) &\
           (x.MS1_Score >= MS1_Score) & (x.vol >= vol) & (x.peptideLens >= peptideLens) & (x.Obs_Mass >= Obs_Mass)


def build_shuffle_seqs(scored_matches, count=20, prefix_len=0, suffix_len=0):
    '''
        :type scored_matches: pd.DataFrame
        :param scored_matches: pd.DataFrame having scored MS2 matches
        :param count: int > 0, number of random peptides to generate per match
        :param prefix_len: int >= 0, number of AA to preserve order on of the
                           start of the match sequence
        :param suffix_len: int >= 0, number of AA to preserve order on of the
                           end of the match sequence
    '''
    split_sequence = scored_matches.Glycopeptide_identifier.str.split(r"([^\[]+)(\[.*\])")
    split_sequence_frame = pd.DataFrame(split_sequence.tolist(), columns=["X", "peptide_sequence", "glycan", "Y"])
    split_sequence_frame["seq_obj"] = split_sequence_frame.peptide_sequence.apply(Sequence)
    shuffles = defaultdict(set)
    for ind, row in split_sequence_frame.iterrows():
        solutions = shuffles[ind]
        mass = scored_matches.Obs_Mass[ind]
        ppm_error = scored_matches.ppm_error[ind]
        ms1_score = scored_matches.MS1_Score[ind]
        glycan = row.glycan
        seq_obj = row.seq_obj
        solutions.add((row.peptide_sequence, glycan, mass, ppm_error, ms1_score))
        while(len(solutions) < count + 1):
            clone = seq_obj.clone()
            pref = seq_obj[:prefix_len]
            suf = seq_obj[-suffix_len:]
            body = seq_obj[prefix_len:len(seq_obj)-suffix_len]
            random.shuffle(body)
            clone.seq = pref + body + suf
            # random.shuffle(clone)
            solutions.add((clone.get_sequence(), glycan, mass, ppm_error, ms1_score))
        solutions.discard((row.peptide_sequence, glycan, mass, ppm_error, ms1_score))
    return shuffles


def calculate_fdr(scored_matches_frame, decoy_matches_frame, predicate_fn=lambda x: x.MS2_Score > 0.5,
                  num_decoys_per_real_mass=20.0):
    pred_passed = scored_matches_frame.apply(predicate_fn, 1)
    decoy_passed = decoy_matches_frame.apply(predicate_fn, 1)
    n_pred_passed = sum(pred_passed)
    n_decoy_passed = sum(decoy_passed)
    total_assignments = n_pred_passed + n_decoy_passed
    fdr = (n_decoy_passed / float(total_assignments)) * (1 + (1 / float(num_decoys_per_real_mass)))
    results_obj = {
        "num_real_matches": n_pred_passed,
        "num_decoy_matches": n_decoy_passed,
        "decoy_to_real_ratio": num_decoys_per_real_mass,
        "false_discovery_rate": fdr,
        "predicted_database": None,
        "decoy_database": None,
        "searched_spectra": None,
        "threshold": (getattr(predicate_fn, "keywords", None))
    }
    return results_obj


def generate_decoy_match_results(scored_matches_path, decon_data, model_file_path,
                                 prefix_len=0, suffix_len=0, ms1_tolerance=1e-5,
                                 ms2_tolerance=2e-5, num_decoys_per_real_mass=20.0,
                                 method="full",
                                 method_init_args=None, method_fit_args=None):
    scored_matches_frame = classify_matches.prepare_model_file(scored_matches_path)
    decoy_file_name = scored_matches_path[:-4] + ".decoy.ion_space.csv"
    random_sequences = build_shuffle_seqs(scored_matches_frame, int(num_decoys_per_real_mass),
                                          prefix_len=prefix_len, suffix_len=suffix_len)
    fh = open(decoy_file_name, "wb")
    fragments = (random_glycopeptide_to_sequence_space(*decoy) for
                 decoy in itertools.chain.from_iterable(random_sequences.values()))
    # Fetch the first fragment set for headers
    frag_sample = fragments.next()
    writer = csv.DictWriter(fh, frag_sample.keys())
    # Write out fragments
    writer.writeheader()
    writer.writerow(frag_sample)
    writer.writerows(fragments)
    fh.close()
    print(decoy_file_name)
    match_file = match_ions.match_frags(decoy_file_name, decon_data, ms1_tolerance, ms2_tolerance)
    print(match_file)
    postprocess_file = postprocess.main(match_file)
    print(postprocess_file)
    classifier = classify_matches.ClassifyTargetWithModelTask(model_file_path, postprocess_file,
                                                              method=method,
                                                              method_init_args=method_init_args,
                                                              method_fit_args=method_fit_args)
    decoy_matches_path = classifier.run()
    return decoy_matches_path


def main(scored_matches_path, decon_data=None, model_file_path=None, decoy_matches_path=None,
         outfile_path=None, num_decoys_per_real_mass=20.0,
         predicate_fns=(make_predicate(MS2_Score=0.5),), prefix_len=0, suffix_len=0,
         ms1_tolerance=1e-5, ms2_tolerance=2e-5, method="full", method_init_args=None,
         method_fit_args=None):

    '''
        Call with deconvolution results and a model to generate decoys and score them, or with
        a pre-existing decoy database.

        :type predicate_fns: Sequence
        :param predicate_fns iterable: containing functions with which to partition both the
        "real" and "decoy" databases. Use `make_predicate` with keyword arguments matching column
        names and numeric thresholds for ease of use and documentation.

        :param outfile_path str: defaults to scored_matches_path[:-4] + "_fdr.csv" will contain the
        resulting FDR statistic for each cutoff.
    '''

    if outfile_path is None:
        outfile_path = scored_matches_path[:-4] + "_fdr.csv"
    if decon_data is not None and model_file_path is not None:
        decoy_matches_path = generate_decoy_match_results(
            scored_matches_path, decon_data, model_file_path, prefix_len=prefix_len,
            suffix_len=suffix_len, ms1_tolerance=ms1_tolerance, ms2_tolerance=ms2_tolerance,
            num_decoys_per_real_mass=num_decoys_per_real_mass,
            method=method, method_init_args=method_init_args, method_fit_args=method_fit_args
        )
    scored_matches_frame = classify_matches.prepare_model_file(scored_matches_path)
    decoy_matches_frame = classify_matches.prepare_model_file(decoy_matches_path)
    results = []
    for predicate in predicate_fns:
        results_obj = calculate_fdr(scored_matches_frame, decoy_matches_frame, predicate, num_decoys_per_real_mass)
        results_obj['predicted_database'] = scored_matches_path
        results_obj['decoy_database'] = decoy_matches_path
        results.append(results_obj)

    results_fh = open(outfile_path, 'wb')
    results_writer = csv.DictWriter(results_fh, results[0].keys())
    results_writer.writeheader()
    results_writer.writerows(results)
    results_fh.close()
    return results


if __name__ == '__main__':
    import argparse
    app = argparse.ArgumentParser()
    app.add_argument("--decon_data", default=None, help="MS2 deconvolution output file path")
    app.add_argument("--model_file", default=None, help="The model to train a classifier on for scoring decoys")
    app.add_argument("-d", "--decoy_matches", default=None, required=False, help="Path to decoy matches\
                     that have already been generated")
    app.add_argument("-s", "--scored_matches", required=True,
                     help="Path to results of matching")
    app.add_argument("-p", "--prefix-length", default=0, required=False, help="Length of peptide prefix to preserve\
            when generating random glycopeptides by shuffling.")
    app.add_argument("-t", "--suffix-length", default=1, required=False, help="Length of peptide suffix to preserve\
            when generating random glycopeptides by shuffling.")
    args = app.parse_args()
    predicates = [make_predicate(MS2_Score=i) for i in [0.2, 0.4, 0.6, 0.8, .9]]
    predicates.extend([make_predicate(MS2_Score=i, peptideLens=10) for i in [0.2, 0.4, 0.6, 0.8, .9]])
    predicates.extend([make_predicate(MS2_Score=i, peptideLens=15) for i in [0.2, 0.4, 0.6, 0.8, .9]])

    main(args.scored_matches, args.decon_data, args.model_file, args.decoy_matches, predicate_fns=predicates)
