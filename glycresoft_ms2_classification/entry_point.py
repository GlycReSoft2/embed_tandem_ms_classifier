import os
import logging
logfile = os.path.expanduser("~/.glycresoft-log")
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()
fh = logging.FileHandler(logfile, mode='w')
fh.setLevel(logging.DEBUG)
logger.addHandler(fh)
logger.addHandler(logging.StreamHandler())

import atexit
import json

# Handles URI Decoding incompatibility with py2K
try:
    import urllib2 as url_parser  # python 2
except ImportError:
    import urllib.parse as url_parser  # python 3

import multiprocessing

from error_code_interface import *
# Import the pipeline modules, wrapping any ImportErrors in the cross-runtime communication
# exception.
try:
    # Pipeline Steps as Modules
    from structure.sequence_space import UnqualifiedModifierException
    from structure.sequence_space import NoSitesFoundException
    from proteomics.msdigest_xml_parser import MSDigestParamters

    import theoretical_glycopeptide
    import match_ions2
    import postprocess2
    from classify_matches import PrepareModelTask
    from classify_matches import ClassifyTargetWithModelTask
    from classify_matches import ModelDiagnosticsTask
    #from classify_matches import CompareModelsDiagnosticTask
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
                                   enzyme_info, n_processes):
    try:
        path = theoretical_glycopeptide.main(ms1_results_file, sites_file, constant_modifications,
                                             variable_modifications, enzyme_info, n_processes)
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
    try:
        task = ClassifyTargetWithModelTask(
            model_file_path, postprocessed_ions_file, method=method, output_path=out,
            method_init_args=method_init_args, method_fit_args=method_fit_args)
        result = task.run()
        return result
    except Exception, e:
        raise ClassificationException(str(e))


def calculate_false_discovery_rate(scored_predictions_file, deconvoluted_spectra, model_file_path,
                                   method="full_random_forest", ms1_match_tolerance=1e-5, ms2_match_tolerance=2e-5,
                                   out=None,  n_decoys=20, predicates=None, n_processes=6):
    try:
        predicates = predicates if predicates is not None else calculate_fdr.default_predicates()
        outfile_path = calculate_fdr.main(scored_matches_path=scored_predictions_file, decon_data=deconvoluted_spectra,
                                          model_file_path=model_file_path, decoy_matches_path=None,
                                          outfile_path=out, num_decoys_per_real_mass=n_decoys,
                                          predicate_fns=predicates, prefix_len=0, suffix_len=0, by_mod_sig=False,
                                          ms1_tolerance=ms1_match_tolerance, ms2_tolerance=ms2_match_tolerance,
                                          method="full_random_forest", method_init_args=None,
                                          method_fit_args=None, n_processes=n_processes)
        return outfile_path
    except NameError, e:
        pass


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
        n_processes=n_processes)
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
    ms1_results_file, glycosylation_sites_file, deconvoluted_spectra_file,
    ms1_match_tolerance, ms2_match_tolerance,
    constant_modification_list=None,
    variable_modification_list=None, enzyme=None,
    method="full_random_forest", model_file=None,
    n_processes=4,
        out=None):
    theoretical_ion_space_file = generate_theoretical_ion_space(
        ms1_results_file, glycosylation_sites_file,
        constant_modification_list, variable_modification_list,
        enzyme, n_processes=n_processes)
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
    classification_results_file = classify_data_by_model(
        postprocessed_ions_file, method=method,
        model_file_path=model_file, out=out)
    print(classification_results_file)
    return classification_results_file


def reclassify_with_model_app_function(target_file, method="full_random_forest", model_file=None, out=None):
    if out is None:
        out = os.path.splitext(target_file)[0] + ".rescored.csv"
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


def main():
    import argparse
    app = argparse.ArgumentParser()
    subparsers = app.add_subparsers()

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

    build_model_app.add_argument("-e", "--enzyme", action="store", help="Name of the enzyme used")
    build_model_app.add_argument("-p", "--protein_prospector_xml", action="store", help="path to msdgist XML file.\
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
                                         action="store", default=None, required=True)
    classify_with_model_app.add_argument(
        "--ms1-match-tolerance", type=float, action="store",
        default=match_ions2.ms1_tolerance_default,
        help="Mass Error Tolerance for matching MS1 masses in PPM")
    classify_with_model_app.add_argument(
        "--ms2-match-tolerance", type=float, action="store",
        default=match_ions2.ms2_tolerance_default,
        help="Mass Error Tolerance for matching MS2 masses in PPM")
    classify_with_model_app.add_argument("--protein-prospector-xml", default=None, help="Parse out modification\
                                 information form the Protein Prospector XML output")
    classify_with_model_app.add_argument(
        "--constant-modification-list", type=str, action="append", default=None,
        help="Pass the list of constant modifications to include in the sequence search space")
    classify_with_model_app.add_argument(
        "--variable-modification-list", type=str, action="append", default=None,
        help="Pass the list of variable modifications to include in the sequence search space")
    classify_with_model_app.add_argument("-e", "--enzyme", action="store", help="Name of the enzyme used")
    classify_with_model_app.add_argument("-p", "--protein_prospector_xml", action="store", help="path to msdgist\
                                XML file. Instead of --enzyme,--constant_modifications and --variable_modifications")
    classify_with_model_app.add_argument("--out", action="store", default=None)
    classify_with_model_app.set_defaults(func=classify_with_model_app_function)

    reclassify_with_model_app = subparsers.add_parser(
        "reclassify-with-model", help="Rerun classification of an matched ion data file")
    reclassify_with_model_app.add_argument("--target-file", action="store", default=None, required=True,
                                           help="Matched ion data file to re-classify")
    reclassify_with_model_app.add_argument(
        "--method", action="store", default="full_random_forest",
        choices=set(ModelDiagnosticsTask.method_table),
        help="Select the model method to use for classification")
    reclassify_with_model_app.add_argument("--model-file",
                                           action="store", default=None, required=True)
    reclassify_with_model_app.set_defaults(func=reclassify_with_model_app_function)

    model_diagnostics_app = subparsers.add_parser(
        "model-diagnostics", help="Given a labeled model, calculate model diagnostics")
    model_diagnostics_app.add_argument(
        "--method", action="store", default="full_random_forest",
        choices=set(ModelDiagnosticsTask.method_table),
        help="Select the model method to use for classification")
    model_diagnostics_app.add_argument("--model-file",
                                       action="store", default=None, required=True)
    model_diagnostics_app.set_defaults(func=model_diagnostics_app_function)
    try:
        # Always generate an error map file on run, for sanity's sake. But users may not have write-access
        # to the install directory. Need to use tempfile or other source?
        # error_map_file = ErrorCodingMeta.build_error_code_map()
        args = app.parse_args()
        args = args.__dict__
        func = args.pop("func")
        debug = args.pop("debug", os.environ.get("GLYCRESOFT_DEBUG", False))

        if 'constant_modification_list' in args:
            args['constant_modification_list'] = uri_decode_list(
                args['constant_modification_list'])

        if 'variable_modification_list' in args:
            args['variable_modification_list'] = uri_decode_list(
                args['variable_modification_list'])

        if 'protein_prospector_xml' in args and args["protein_prospector_xml"] is not None:
            ms_digest = MSDigestParamters.parse(args["protein_prospector_xml"])
            args["constant_modification_list"] = ms_digest.constant_modifications
            args["variable_modification_list"] = ms_digest.variable_modifications
            args["enzyme"] = ms_digest.enzyme
        args.pop("protein_prospector_xml", None)

        param_file = args.pop("parameter_file", None)
        if param_file is not None:
            params = load_parameters_from_json(param_file)
            args.update(params)
        if debug:
            logger.info(json.dumps(args, indent=4))
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
