import os
import sys

# Ensure our scripts are imported before any others since we can't depend upon
# the packaging middleware to handle relative imports, or should we vendorize
# pip or easy_install?
sys.path = [(os.path.dirname(os.path.abspath(__file__)))] + sys.path


#import atexit
import json

# Handles URI Decoding incompatibility with py2K
try:
    import urllib2 as url_parser  # python 2
except ImportError:
    import urllib.parse as url_parser  # python 3

# Import all error code exceptions. Yes. All of them. A valid use of *
from error_code_interface import *
# Import the pipeline modules, wrapping any ImportErrors in the cross-runtime communication
# exception.
try:
    # Pipeline Steps as Modules
    from structure.sequence_space import UnqualifiedModifierException
    from structure.sequence_space import NoSitesFoundException

    from protein_prospector.xml_parser import MSDigestParamters

    import theoretical_glycopeptide
    import match_ions
    import postprocess
    from classify_matches import PrepareModelTask
    from classify_matches import ClassifyTargetWithModelTask
    from classify_matches import ModelDiagnosticsTask
except ImportError, e:
    exception = ImportErrorWrapperException(str(e))
    print(exception)
    exit(exception.errcode)


def load_parameters(param_file):
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


def generate_theoretical_ion_space(ms1_results_file, sites_file, constant_modifications, variable_modifications):
    try:
        return theoretical_glycopeptide.main(ms1_results_file, sites_file,
                                             constant_modifications, variable_modifications, None)
    except NoSitesFoundException, e:
        raise NoSitesFoundWrapperException(str(e))
    except UnqualifiedModifierException, e:
        raise UnqualifiedModificationWrapperException(str(e))


def match_deconvoluted_ions(theoretical_ion_space, deconvoluted_spectra, ms1_match_tolerance, ms2_match_tolerance):
    try:
        return match_ions.match_frags(
            theoretical_ion_space, deconvoluted_spectra,
            ms1_match_tolerance, ms2_match_tolerance, split_decon_data=True,
            outfile=None)
    except match_ions.NoIonsMatchedException, e:
        raise NoIonsMatchedException(str(e))


def postprocess_matches(matched_ions_file):
    return postprocess.main(matched_ions_file)


def prepare_model_file(postprocessed_ions_file, method="full", out=None):
    try:
        task = PrepareModelTask(
            model_file_path=postprocessed_ions_file, method="full", output_path=out)
        result = task.run()
        return result
    except Exception, e:
        raise ModelFitException(str(e))


def classify_data_by_model(
    postprocessed_ions_file, model_file_path, method="full", out=None, method_init_args=None,
        method_fit_args=None):
    try:
        task = ClassifyTargetWithModelTask(
            model_file_path, postprocessed_ions_file, method=method, output_path=out,
            method_init_args=method_init_args, method_fit_args=method_fit_args)
        result = task.run()
        return result
    except Exception, e:
        raise ClassificationException(str(e))


def windows_to_unix_path(path):
    return path.replace("\\", "/")


def clean_up_files(*files):
    for f in files:
        try:
            os.remove(f)
        except:
            pass

intermediary_files = []


def __intermediary_files():
    print(intermediary_files)


def uri_decode(uri):
    return url_parser.unquote(uri)


def uri_decode_list(uri_list=None):
    print(uri_list)
    if uri_list is None:
        return None
    return map(uri_decode, uri_list)


def build_model_app_function(
    ms1_results_file, glycosylation_sites_file, deconvoluted_spectra_file,
    ms1_match_tolerance, ms2_match_tolerance,
    constant_modification_list=None,
    variable_modification_list=None,
    method="full",
        out=None):
    theoretical_ion_space_file = generate_theoretical_ion_space(
        ms1_results_file, glycosylation_sites_file,
        constant_modification_list, variable_modification_list)
    # print(theoretical_ion_space_file)
    intermediary_files.append(theoretical_ion_space_file)
    matched_ions_file = match_deconvoluted_ions(
        theoretical_ion_space_file, deconvoluted_spectra_file, ms1_match_tolerance, ms2_match_tolerance)
    # print(matched_ions_file)
    intermediary_files.append(matched_ions_file)
    postprocessed_ions_file = postprocess_matches(matched_ions_file)
    # print(postprocessed_ions_file)
    intermediary_files.append(postprocessed_ions_file)
    model_file = prepare_model_file(
        postprocessed_ions_file, method=method, out=out)
    print(model_file)


def classify_with_model_app_function(
    ms1_results_file, glycosylation_sites_file, deconvoluted_spectra_file,
    ms1_match_tolerance, ms2_match_tolerance,
    constant_modification_list=None,
    variable_modification_list=None,
    method="full", model_file=None,
        out=None):
    theoretical_ion_space_file = generate_theoretical_ion_space(
        ms1_results_file, glycosylation_sites_file,
        constant_modification_list, variable_modification_list)
    # print(theoretical_ion_space_file)
    intermediary_files.append(theoretical_ion_space_file)
    matched_ions_file = match_deconvoluted_ions(
        theoretical_ion_space_file, deconvoluted_spectra_file, ms1_match_tolerance, ms2_match_tolerance)
    # print(matched_ions_file)
    intermediary_files.append(matched_ions_file)
    postprocessed_ions_file = postprocess_matches(matched_ions_file)
    # print(postprocessed_ions_file)
    intermediary_files.append(postprocessed_ions_file)
    classification_results_file = classify_data_by_model(
        postprocessed_ions_file, method=method,
        model_file_path=model_file, out=out)
    print(classification_results_file)


def model_diagnostics_app_function(model_file, method="full"):
    task = ModelDiagnosticsTask(model_file_path=model_file, method=method)
    result = task.run()
    print(result)


def main():
    #atexit.register(lambda: clean_up_files(*intermediary_files))
    import argparse
    app = argparse.ArgumentParser()
    subparsers = app.add_subparsers()

    app.add_argument(
        "-d", "--debug", action='store_true', default=False, required=False)

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
        "--deconvoluted-spectra-file", action="store", required=True)
    build_model_app.add_argument("--method", action="store", default="full",
                                 choices=set(PrepareModelTask.method_table),
                                 help="Select the model method to use for classification")
    build_model_app.add_argument(
        "--ms1-match-tolerance", type=float, action="store",
        default=match_ions.ms1_tolerance_default,
        help="Mass Error Tolerance for matching MS1 masses in PPM")
    build_model_app.add_argument(
        "--ms2-match-tolerance", type=float, action="store",
        default=match_ions.ms2_tolerance_default,
        help="Mass Error Tolerance for matching MS2 masses in PPM")
    build_model_app.add_argument("--protein-prospector-xml", default=None, help="Parse out modification\
                                 information form the Protein Prospector XML output")
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
        "--method", action="store", default="full",
        choices=set(ClassifyTargetWithModelTask.method_table),
        help="Select the model method to use for classification")
    classify_with_model_app.add_argument("--model-file",
                                         action="store", default=None, required=True)
    classify_with_model_app.add_argument(
        "--ms1-match-tolerance", type=float, action="store",
        default=match_ions.ms1_tolerance_default,
        help="Mass Error Tolerance for matching MS1 masses in PPM")
    classify_with_model_app.add_argument(
        "--ms2-match-tolerance", type=float, action="store",
        default=match_ions.ms2_tolerance_default,
        help="Mass Error Tolerance for matching MS2 masses in PPM")
    classify_with_model_app.add_argument("--protein-prospector-xml", default=None, help="Parse out modification\
                                 information form the Protein Prospector XML output")
    classify_with_model_app.add_argument(
        "--constant-modification-list", type=str, action="append", default=None,
        help="Pass the list of constant modifications to include in the sequence search space")
    classify_with_model_app.add_argument(
        "--variable-modification-list", type=str, action="append", default=None,
        help="Pass the list of variable modifications to include in the sequence search space")

    classify_with_model_app.add_argument("--out", action="store", default=None)
    classify_with_model_app.set_defaults(func=classify_with_model_app_function)

    model_diagnostics_app = subparsers.add_parser(
        "model-diagnostics", help="Given a labeled model, calculate model diagnostics")
    model_diagnostics_app.add_argument(
        "--method", action="store", default="full",
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
        debug = args.pop("debug", False)

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
        args.pop("protein_prospector_xml", None)

        param_file = args.pop("parameter_file", None)
        if param_file is not None:
            params = load_parameters(param_file)
            args.update(params)
        if debug:
            print(json.dumps(args, indent=4))
        func(**args)
    except GlycReSoftInterprocessCommunicationException, e:
        print(e)
        exit(e.errcode)
    exit(0)

if __name__ == '__main__':
    main()
