import logging
logger = logging.getLogger("calculate_fdr")

from glycresoft_ms2_classification import classify_matches
from glycresoft_ms2_classification import match_ions2
from glycresoft_ms2_classification import postprocess2
from glycresoft_ms2_classification.classify_matches import ClassifyTargetWithModelTask
from glycresoft_ms2_classification.prediction_tools.false_discovery_rate import make_decoys
from glycresoft_ms2_classification.prediction_tools.false_discovery_rate import parameter_search_optimizer as optimize_fdr


def predicate_base(x, MS2_Score=0, meanCoverage=0, meanHexNAcCoverage=0, percentUncovered=1, MS1_Score=0,
                   peptideLens=0, Obs_Mass=0, numStubs=-1):
    return (x.MS2_Score >= MS2_Score) & (x.meanCoverage >= meanCoverage) & (x.percentUncovered <= percentUncovered) &\
           (x.MS1_Score >= MS1_Score) & (x.peptideLens >= peptideLens) & (x.Obs_Mass >= Obs_Mass) &\
           (x.numStubs >= numStubs) & (x.meanHexNAcCoverage >= meanHexNAcCoverage)


def calculate_fdr(scored_matches_frame, decoy_matches_frame, predicate_fn=lambda x: x.MS2_Score > 0.5,
                  num_decoys_per_real_mass=20.0):
    pred_passed = scored_matches_frame.apply(predicate_fn, 1)
    decoy_passed = decoy_matches_frame.apply(predicate_fn, 1)
    n_pred_passed = sum(pred_passed)
    if len(decoy_matches_frame.index) == 0:
        n_decoy_passed = -1
        fdr = 0
    else:
        n_decoy_passed = sum(decoy_passed)
        total_assignments = n_pred_passed + n_decoy_passed
        fdr = (n_decoy_passed / float(total_assignments)) * \
            (1 + (1 / float(num_decoys_per_real_mass)))

    if n_pred_passed == 0 and n_decoy_passed <= 0:
        return None

    results_obj = {
        "num_real_matches": n_pred_passed,
        "num_decoy_matches": n_decoy_passed,
        "decoy_to_real_ratio": num_decoys_per_real_mass,
        "false_discovery_rate": fdr,
        "threshold": (getattr(predicate_fn, "sig", None))
    }
    return results_obj


def generate_decoy_match_results(scored_matches_path, decon_data, model_file_path,
                                 prefix_len=0, suffix_len=0, ms1_tolerance=1e-5,
                                 ms2_tolerance=2e-5, num_decoys_per_real_mass=20.0,
                                 method="full_random_forest", random_only=False,
                                 method_init_args=None, method_fit_args=None,
                                 n_processes=6, outfile_path=None):

    decoy_file_name = make_decoys.taskmain(scored_matches_path,
                                           count=int(num_decoys_per_real_mass),
                                           prefix_len=prefix_len, suffix_len=suffix_len,
                                           n_processes=n_processes, random_only=random_only, out=outfile_path)
    logger.info("Decoy Ion Space: %s", decoy_file_name)
    match_ions2.match_frags(
        decoy_file_name, decon_data,
        ms1_tolerance, ms2_tolerance, n_processes=n_processes)
    logger.info("Decoy Matches Done")
    postprocess_file, postprocess_data = postprocess2.main(decoy_file_name)
    logger.info("Decoys Postprocessed: %s", postprocess_file)
    classifier = classify_matches.ClassifyTargetWithModelTask(model_file_path, postprocess_file,
                                                              method=method,
                                                              method_init_args=method_init_args,
                                                              method_fit_args=method_fit_args)
    decoy_matches_path = classifier.run()
    return decoy_matches_path


def main(scored_matches_path, decon_data=None, model_file_path=None, decoy_matches_path=None,
         outfile_path=None, num_decoys_per_real_mass=20.0, random_only=False,
         predicate_fns=None, prefix_len=0, suffix_len=0, by_mod_sig=False,
         ms1_tolerance=None, ms2_tolerance=None, method="full_random_forest", method_init_args=None,
         method_fit_args=None, n_processes=6):
    '''
        Call with deconvolution results and a model to generate decoys and score them, or with
        a pre-existing decoy database.

        :type predicate_fns: Sequence
        :param predicate_fns iterable: containing functions with which to partition both the
        "real" and "decoy" databases. Use `make_predicate` with keyword arguments matching column
        names and numeric thresholds for ease of use and documentation.

        :param outfile_path str: defaults to scored_matches_path[:-4] + "_fdr.json" will contain the
        resulting FDR statistic for each cutoff.
    '''
    scored_matches_frame = classify_matches.prepare_model_file(scored_matches_path)
    decoy_matches_frame = None
    if ms1_tolerance is None:
        ms1_tolerance = scored_matches_frame.metadata["ms1_ppm_tolerance"]
    if ms2_tolerance is None:
        ms2_tolerance = scored_matches_frame.metadata["ms2_ppm_tolerance"]

    if outfile_path is None:
        outfile_path = scored_matches_path[:-5] + "_fdr.json"
    if decon_data is not None and model_file_path is not None:
        decoy_matches_path = generate_decoy_match_results(
            scored_matches_path, decon_data, model_file_path, prefix_len=prefix_len,
            suffix_len=suffix_len, ms1_tolerance=ms1_tolerance, ms2_tolerance=ms2_tolerance,
            num_decoys_per_real_mass=num_decoys_per_real_mass, random_only=random_only,
            method=method, method_init_args=method_init_args, method_fit_args=method_fit_args,
            n_processes=n_processes, outfile_path=outfile_path)
        decoy_matches_frame = classify_matches.prepare_model_file(
            decoy_matches_path)
    elif model_file_path is not None and decoy_matches_path is not None:
        scored_matches_frame = ClassifyTargetWithModelTask(
            model_file_path, scored_matches_path, method=method).run(False)
        decoy_matches_frame = ClassifyTargetWithModelTask(
            model_file_path, decoy_matches_path, method=method).run(False)
    else:
        scored_matches_frame = classify_matches.prepare_model_file(
            scored_matches_path)
        decoy_matches_frame = classify_matches.prepare_model_file(
            decoy_matches_path)
    logger.info("Evaluating predicates")
    fdr_search = optimize_fdr.CountExclusion(
        scored_matches_frame, decoy_matches_frame,
        decoy_matches_frame.metadata["decoy_ratio"], ["MS2_Score"])
    fdr_search.optimize()

    scored_matches_frame.metadata["fdr"] = fdr_search.compress()
    scored_matches_frame["call"] = scored_matches_frame.Glycopeptide_identifier.isin(
        scored_matches_frame.kvquery().Glycopeptide_identifier)
    logger.info("Accepted %d predictions with FDR %f",
                scored_matches_frame.call.sum(),
                scored_matches_frame.optimize_fdr().iloc[0]["false_discovery_rate"])
    scored_matches_frame.serialize(outfile_path)

    return outfile_path


def default_predicates():
    predicates = None
    return predicates


if __name__ == '__main__':
    import argparse
    app = argparse.ArgumentParser()
    app.add_argument("--decon_data", default=None, help="MS2 deconvolution output file path")
    app.add_argument("--model_file", default=None, help="The model to train a classifier on for\
     scoring decoys")
    app.add_argument("-d", "--decoy_matches", default=None, required=False, help="Path to decoy matches\
                     that have already been generated")
    app.add_argument("-s", "--scored_matches", required=True,
                     help="Path to results of matching")
    app.add_argument("-p", "--prefix-length", default=0, required=False, help="Length of peptide prefix to preserve\
            when generating random glycopeptides by shuffling.")
    app.add_argument("-t", "--suffix-length", default=1, required=False, help="Length of peptide suffix to preserve\
            when generating random glycopeptides by shuffling.")
    app.add_argument("--modification_signature", required=False, default=False,
                     action="store_true", help="Should FDR by computed by modification signature groups?")
    args = app.parse_args()
    predicates = default_predicates()

    main(args.scored_matches, args.decon_data, args.model_file,
         args.decoy_matches, predicate_fns=predicates, by_mod_sig=args.modification_signature)
