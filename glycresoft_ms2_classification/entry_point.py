import os
# import atexit
import json
import logging
try:
    import urllib2 as url_parser  # python 2
except ImportError:
    import urllib.parse as url_parser  # python 3
import multiprocessing
from error_code_interface import (NoSitesFoundWrapperException, ImportErrorWrapperException,
                                  UnqualifiedModificationWrapperException, NoIonsMatchedException,
                                  ModelFitException)
try:
    logfile = "./.glycresoft-log"
    open(logfile, 'w').close()
except Exception, e:
    print(e)
    logfile = os.path.expanduser("~/.glycresoft-log")
    print("logging to ~/")

logging.basicConfig(level=logging.DEBUG, filename=logfile, filemode='w',
                    format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                    datefmt="%H:%M:%S")
logger = logging.getLogger()

from error_code_interface import *

# Handles URI Decoding incompatibility with py2K
# Import the pipeline modules, wrapping any ImportErrors in the cross-runtime communication
# exception.
try:
    from glycresoft_ms2_classification.utils import config_loader
    # Pipeline Steps as Modules
    from glycresoft_ms2_classification.structure.sequence_space import UnqualifiedModifierException
    from glycresoft_ms2_classification.structure.sequence_space import NoSitesFoundException
    from glycresoft_ms2_classification.proteomics.msdigest_xml_parser import MSDigestParamters
    from glycresoft_ms2_classification import theoretical_glycopeptide
    from glycresoft_ms2_classification import match_ions2
    from glycresoft_ms2_classification import postprocess2
    from glycresoft_ms2_classification import calculate_fdr
    from glycresoft_ms2_classification.classify_matches import PrepareModelTask
    from glycresoft_ms2_classification.classify_matches import ClassifyTargetWithModelTask
    from glycresoft_ms2_classification.classify_matches import ModelDiagnosticsTask
except ImportError, e:
    exception = ImportErrorWrapperException(str(e))
    print(exception)
    exit(exception.errcode)


def load_parameters_from_json(param_file):
    handle = None
    if hasattr(param_file, 'read'):
        handle = param_file
    elif isinstance(param_file, str):
        try:
            handle = open(param_file, 'rb')
        except:
            raise IOError("Could not read from parameter file")
    params = json.load(handle)
    return params


def generate_theoretical_ion_space(ms1_results_file, sites_file, constant_modifications, variable_modifications,
                                   enzyme_info, n_processes, out=None):
    try:
        path = theoretical_glycopeptide.main(ms1_results_file, sites_file, constant_modifications,
                                             variable_modifications, enzyme_info, n_processes, output_file=out)
        return path
    except NoSitesFoundException, e:
        raise NoSitesFoundWrapperException(str(e))
    except UnqualifiedModifierException, e:
        raise UnqualifiedModificationWrapperException(str(e))


def match_deconvoluted_ions(theoretical_ion_space, deconvoluted_spectra,
                            ms1_match_tolerance, ms2_match_tolerance, n_processes):
    try:
        path = match_ions2.match_frags(
            theoretical_ion_space, deconvoluted_spectra,
            ms1_match_tolerance, ms2_match_tolerance, split_decon_data=True,
            n_processes=n_processes,
            outfile=None)
        return path
    except match_ions2.NoIonsMatchedException, e:
        raise NoIonsMatchedException(str(e))


def postprocess_matches(matched_ions_file):
    path, data = postprocess2.main(matched_ions_file)
    return path


def prepare_model_file(postprocessed_ions_file, method="full_random_forest", out=None):
    try:
        task = PrepareModelTask(
            model_file_path=postprocessed_ions_file, method="full_random_forest", output_path=out)
        result = task.run()
        return result
    except Exception, e:
        raise ModelFitException(str(e))


def classify_data_by_model(
    postprocessed_ions_file, model_file_path, method="full_random_forest", out=None, method_init_args=None,
        method_fit_args=None):
    task = ClassifyTargetWithModelTask(
        model_file_path, postprocessed_ions_file, method=method, output_path=out,
        method_init_args=method_init_args, method_fit_args=method_fit_args)
    result = task.run()
    return result


def calculate_false_discovery_rate(scored_predictions_file, deconvoluted_spectra=None, model_file_path=None,
                                   decoy_matches_path=None,
                                   method="full_random_forest", ms1_match_tolerance=1e-5, ms2_match_tolerance=2e-5,
                                   out=None,  n_decoys=1, predicates=None, random_only=False,
                                   n_processes=6, prefix_len=0, suffix_len=1):
    predicates = predicates if predicates is not None else calculate_fdr.default_predicates()
    outfile_path = calculate_fdr.main(scored_matches_path=scored_predictions_file, decon_data=deconvoluted_spectra,
                                      model_file_path=model_file_path, decoy_matches_path=decoy_matches_path,
                                      outfile_path=out, num_decoys_per_real_mass=n_decoys,
                                      predicate_fns=predicates, prefix_len=prefix_len, suffix_len=suffix_len,
                                      by_mod_sig=False, random_only=random_only,
                                      ms1_tolerance=ms1_match_tolerance, ms2_tolerance=ms2_match_tolerance,
                                      method=method, method_init_args=None,
                                      method_fit_args=None, n_processes=n_processes)
    return outfile_path


def build_model_app_function(
    ms1_results_file, glycosylation_sites_file, deconvoluted_spectra_file,
    ms1_match_tolerance, ms2_match_tolerance,
    constant_modification_list=None,
    variable_modification_list=None, enzyme=None,
    method="full_random_forest",
    n_processes=4,
        out=None):
    theoretical_ion_space_file = generate_theoretical_ion_space(
        ms1_results_file, glycosylation_sites_file,
        constant_modification_list, variable_modification_list, enzyme,
        out=out, n_processes=n_processes)
    # print(theoretical_ion_space_file)
    intermediary_files.append(theoretical_ion_space_file)
    matched_ions_file = match_deconvoluted_ions(
        theoretical_ion_space_file, deconvoluted_spectra_file,
        ms1_match_tolerance, ms2_match_tolerance,
        n_processes=n_processes)
    # print(matched_ions_file)
    intermediary_files.append(matched_ions_file)
    postprocessed_ions_file = postprocess_matches(matched_ions_file)
    # print(postprocessed_ions_file)
    intermediary_files.append(postprocessed_ions_file)
    model_file = prepare_model_file(
        postprocessed_ions_file, method=method, out=out)
    print(model_file)
    return model_file


def classify_with_model_app_function(
    ms1_results_file=None, glycosylation_sites_file=None, deconvoluted_spectra_file=None,
    ms1_match_tolerance=None, ms2_match_tolerance=None,
    constant_modification_list=None, prefix_length=0, suffix_length=1,
    variable_modification_list=None, enzyme=None,
    method="full_random_forest", model_file=None,
    n_processes=4, random_only=False, decoy_to_real_ratio=1,
    out=None, search_space_db_path=None):
    if search_space_db_path is None:
        theoretical_ion_space_file = generate_theoretical_ion_space(
            ms1_results_file, glycosylation_sites_file,
            constant_modification_list, variable_modification_list,
            enzyme, out=out, n_processes=n_processes)
        intermediary_files.append(theoretical_ion_space_file)
    else:
        theoretical_ion_space_file = search_space_db_path
    matched_ions_file = match_deconvoluted_ions(
        theoretical_ion_space_file, deconvoluted_spectra_file,
        ms1_match_tolerance, ms2_match_tolerance,
        n_processes=n_processes)
    intermediary_files.append(matched_ions_file)
    postprocessed_ions_file = postprocess_matches(matched_ions_file)
    intermediary_files.append(postprocessed_ions_file)
    classification_results_file = classify_data_by_model(
        postprocessed_ions_file, method=method,
        model_file_path=model_file, out=out)
    logger.info("Calculating FDR")
    fdr_results_file = calculate_false_discovery_rate(
        classification_results_file, deconvoluted_spectra=deconvoluted_spectra_file,
        model_file_path=model_file,
        method="full_random_forest", random_only=random_only, prefix_len=prefix_length,
        suffix_len=suffix_length,
        ms1_match_tolerance=ms1_match_tolerance, ms2_match_tolerance=ms2_match_tolerance,
        out=None,  n_decoys=decoy_to_real_ratio, n_processes=n_processes)
    print(fdr_results_file)
    return fdr_results_file


def reclassify_with_model_app_function(target_file, method="full_random_forest",
                                       model_file=None, out=None, n_processes=6):
    if out is None:
        out = os.path.splitext(target_file)[0] + ".rescored"
    classification_results_file = classify_data_by_model(
        target_file, method=method,
        model_file_path=model_file, out=out)
    print(classification_results_file)
    return classification_results_file


def model_diagnostics_app_function(model_file, method="full_random_forest", **kwargs):
    task = ModelDiagnosticsTask(model_file_path=model_file, method=method)
    result = task.run()
    print(result)
    return result


def calculate_fdr_app_function(predictions_file, model_file=None, decoys_file=None,
                               deconvoluted_spectra_file=None, prefix_length=0, suffix_length=0,
                               decoy_to_real_ratio=1, method="full_random_forest",
                               predicate_fns=None, random_only=False,
                               out=None, n_processes=6):
    if (decoys_file is None and deconvoluted_spectra_file is None) or\
       (deconvoluted_spectra_file is not None and model_file is None):
        print("You must either specify a file containing pre-computed and matched decoys (--decoy-file), \
or a deconvoluted tandem spectra file and model to match and score newly generated decoys with\
 (--deconvoluted-spectra-file, --model-file).")

    fdr_result = calculate_false_discovery_rate(predictions_file, model_file_path=model_file,
                                                decoy_matches_path=decoys_file,
                                                deconvoluted_spectra=deconvoluted_spectra_file,
                                                method=method, out=out, n_decoys=decoy_to_real_ratio,
                                                random_only=random_only, prefix_len=prefix_length,
                                                suffix_len=suffix_length,
                                                predicates=None, n_processes=n_processes)
    print(fdr_result)
    return fdr_result


def clean_up_files(*files):
    for f in files:
        try:
            os.remove(f)
        except:
            pass

intermediary_files = []


def __intermediary_files():
    os.remove(*intermediary_files)


def uri_decode(uri):
    return url_parser.unquote(uri)


def uri_decode_list(uri_list=None):
    if uri_list is None:
        return None
    return map(uri_decode, uri_list)


def main():
    import argparse
    app = argparse.ArgumentParser()
    subparsers = app.add_subparsers()
    app.add_argument("-c", "--config", type=str, default=None, required=False, help="Path to configuration file")
    app.add_argument(
        "-d", "--debug", action='store_true', default=False, required=False)
    app.add_argument("-n", "--n-processes", type=int, action="store",
                     default=4, help="Number of procresses to use")

    # BUILD MODEL
    build_model_app = subparsers.add_parser(
        "build-model", help="Build a model and prepare it to be labeled")
    build_model_app.add_argument(
        "--parameter-file", action="store", default=None)
    build_model_app.add_argument(
        "--ms1-results-file", action="store", required=True)
    build_model_app.add_argument(
        "--glycosylation-sites-file", action="store", required=True)

    build_model_app.add_argument(
        "-e", "--enzyme", action="store", help="Name of the enzyme used")
    build_model_app.add_argument("-p", "--protein-prospector-xml", action="store", help="path to msdgist XML file.\
     Instead of --enzyme,--constant_modifications and --variable_modifications")

    build_model_app.add_argument(
        "--deconvoluted-spectra-file", action="store", required=True)

    build_model_app.add_argument("--method", action="store", default="full_random_forest",
                                 choices=set(PrepareModelTask.method_table),
                                 help="Select the model method to use for classification")
    build_model_app.add_argument(
        "--ms1-match-tolerance", type=float, action="store",
        default=match_ions2.ms1_tolerance_default,
        help="Mass Error Tolerance for matching MS1 masses in PPM")

    build_model_app.add_argument(
        "--ms2-match-tolerance", type=float, action="store",
        default=match_ions2.ms2_tolerance_default,
        help="Mass Error Tolerance for matching MS2 masses in PPM")

    build_model_app.add_argument(
        "--constant-modification-list", type=str, action="append", default=None,
        help="Pass the list of constant modifications to include in the sequence search space")
    build_model_app.add_argument(
        "--variable-modification-list", type=str, action="append", default=None,
        help="Pass the list of variable modifications to include in the sequence search space")

    build_model_app.add_argument("--out", action="store", default=None)
    build_model_app.set_defaults(func=build_model_app_function)

    # CLASSIFY WITH MODEL
    classify_with_model_app = subparsers.add_parser(
        "classify-with-model", help="Classify a data set using a labeled model")
    classify_with_model_app.add_argument(
        "--parameter-file", action="store", default=None)
    classify_with_model_app.add_argument(
        "--ms1-results-file", action="store", required=True)
    classify_with_model_app.add_argument(
        "--glycosylation-sites-file", action="store", required=True)
    classify_with_model_app.add_argument(
        "--deconvoluted-spectra-file", action="store", required=True)
    classify_with_model_app.add_argument(
        "--method", action="store", default="full_random_forest",
        choices=set(ClassifyTargetWithModelTask.method_table),
        help="Select the model method to use for classification")
    classify_with_model_app.add_argument("--model-file",
                                         action="store", default="naive", required=True)
    classify_with_model_app.add_argument(
        "--ms1-match-tolerance", type=float, action="store",
        default=match_ions2.ms1_tolerance_default,
        help="Mass Error Tolerance for matching MS1 masses in PPM")
    classify_with_model_app.add_argument(
        "--ms2-match-tolerance", type=float, action="store",
        default=match_ions2.ms2_tolerance_default,
        help="Mass Error Tolerance for matching MS2 masses in PPM")
    classify_with_model_app.add_argument(
        "--constant-modification-list", type=str, action="append", default=None,
        help="Pass the list of constant modifications to include in the sequence search space")
    classify_with_model_app.add_argument(
        "--variable-modification-list", type=str, action="append", default=None,
        help="Pass the list of variable modifications to include in the sequence search space")
    classify_with_model_app.add_argument(
        "-e", "--enzyme", action="store", help="Name of the enzyme used")
    classify_with_model_app.add_argument("-p", "--protein-prospector-xml", action="store", help="path to msdgist\
                                XML file. Instead of --enzyme,--constant_modifications and --variable_modifications")
    classify_with_model_app.add_argument("--out", action="store", default=None)
    classify_with_model_app.set_defaults(func=classify_with_model_app_function)

    classify_with_model_app.add_argument("--decoy-to-real-ratio", action="store", default=1, type=int, help="Number of\
        decoys per prediction sequence")
    classify_with_model_app.add_argument("--random-only", action="store_true", default=False, help="Don't\
        generate shuffled decoys, only randomized sequences")
    classify_with_model_app.add_argument("--prefix-length", default=0, required=False, type=int,
                                         help="Length of peptide prefix to preserve when generating\
                                          random glycopeptides by shuffling.")
    classify_with_model_app.add_argument("--suffix-length", default=1, required=False, type=int,
                                         help="Length of peptide suffix to preserve when generating\
                                          random glycopeptides by shuffling.")

    reclassify_with_model_app = subparsers.add_parser(
        "reclassify-with-model", help="Rerun classification of an matched ion data file")
    reclassify_with_model_app.add_argument("--target-file", action="store", default=None, required=True,
                                           help="Matched ion data file to re-classify")
    reclassify_with_model_app.add_argument(
        "--method", action="store", default="full_random_forest",
        choices=set(ModelDiagnosticsTask.method_table),
        help="Select the model method to use for classification")
    reclassify_with_model_app.add_argument("--model-file",
                                           action="store", default="naive", required=True)
    reclassify_with_model_app.set_defaults(
        func=reclassify_with_model_app_function)

    # Simple Diagnostic plots for testing a model on itself
    model_diagnostics_app = subparsers.add_parser(
        "model-diagnostics", help="Given a labeled model, calculate model diagnostics")
    model_diagnostics_app.add_argument(
        "--method", action="store", default="full_random_forest",
        choices=set(ModelDiagnosticsTask.method_table),
        help="Select the model method to use for classification")
    model_diagnostics_app.add_argument("--model-file",
                                       action="store", default="naive", required=True)
    model_diagnostics_app.set_defaults(func=model_diagnostics_app_function)

    # Stand alone False Discovery Rate Calculations. Either reuse decoys or create them anew. May also
    # include re-scoring for the predictions and decoys
    calculate_fdr_app = subparsers.add_parser("calculate-fdr", help="Given a set of predictions from a collection\
     of data, estimate the false discovery rate")
    calculate_fdr_app.add_argument("--predictions-file", required=True, help="Path to predictions file generated by\
     classify-with-model, build-model, or reclassify-with-model")
    calculate_fdr_app.add_argument(
        "--deconvoluted-spectra-file", action="store", default=None)
    calculate_fdr_app.add_argument(
        "--method", action="store", default="full_random_forest",
        choices=set(ModelDiagnosticsTask.method_table),
        help="Select the model method to use for classification")
    calculate_fdr_app.add_argument(
        "--decoys-file", default=None, help="A file containing precomputed decoy sequence matches")
    calculate_fdr_app.add_argument("--decoy-to-real-ratio", action="store", default=1, type=int, help="Number of\
        decoys per prediction sequence")
    calculate_fdr_app.add_argument("--random-only", action="store_true", default=False, help="Don't\
        generate shuffled decoys, only randomized sequences")
    calculate_fdr_app.add_argument("--prefix-length", default=0, required=False, type=int,
                                   help="Length of peptide prefix to preserve when generating\
                                          random glycopeptides by shuffling.")
    calculate_fdr_app.add_argument("--suffix-length", default=1, required=False, type=int,
                                   help="Length of peptide suffix to preserve when generating\
                                          random glycopeptides by shuffling.")
    calculate_fdr_app.add_argument("--model-file",
                                   action="store", default="naive", required=False)
    calculate_fdr_app.add_argument("--out", action="store", default=None)
    calculate_fdr_app.set_defaults(func=calculate_fdr_app_function)

    try:
        args = app.parse_args()
        args = args.__dict__
        func = args.pop("func")
        debug = args.pop("debug", os.environ.get("GLYCRESOFT_DEBUG", False))
        config_path = args.pop("config")
        if config_path is not None:
            config_loader.load(config_path)
        logger.debug("Config: %r", json.dumps(config_loader.gather(), indent=4))
        if 'constant_modification_list' in args:
            args['constant_modification_list'] = uri_decode_list(
                args['constant_modification_list'])

        if 'variable_modification_list' in args:
            args['variable_modification_list'] = uri_decode_list(
                args['variable_modification_list'])

        if 'protein_prospector_xml' in args and args["protein_prospector_xml"] is not None:
            ms_digest = MSDigestParamters.parse(args["protein_prospector_xml"])
            args[
                "constant_modification_list"] = ms_digest.constant_modifications
            args[
                "variable_modification_list"] = ms_digest.variable_modifications
            args["enzyme"] = ms_digest.enzyme
        args.pop("protein_prospector_xml", None)

        param_file = args.pop("parameter_file", None)
        if param_file is not None:
            params = load_parameters_from_json(param_file)
            args.update(params)

        if args['out'] is not None and args['out'][0] != os.sep:
            args['out'] = ".{0}{1}".format(os.sep, args['out'])

        logger.info(json.dumps(args, indent=4))
        if debug:
            pass
        else:
            #atexit.register(lambda: clean_up_files(*intermediary_files))
            pass
        logger.debug("Entering main program")
        func(**args)
    except GlycReSoftInterprocessCommunicationException, e:
        logger.debug("An error occurred", exc_info=e)
        exit(e.errcode)
    except MemoryError, e:
        logger.debug("An error occurred", exc_info=e)
        exit(MemoryErrorWrapperException.errcode)
    except Exception, e:
        logger.debug("An error occurred", exc_info=e)
        exit(-255)
    exit(0)

if __name__ == '__main__':
    multiprocessing.freeze_support()
    main()
