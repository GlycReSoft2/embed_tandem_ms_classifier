import atexit
import json
import os
import sys
import subprocess

# Handles URI Decoding incompatibility with py2K
try:
    import urllib2 as url_parser  # python 2
except ImportError:
    import urllib.parse as url_parser  # python 3

# Ensure our scripts are imported before any others since we can't depend upon
# the packaging middleware to handle relative imports, or should we vendorize
# pip or easy_install?
sys.path = [(os.path.dirname(os.path.abspath(__file__)))] + sys.path

# Pipeline Steps as Modules
import theoretical_glycopeptide
import match_ions
import postprocess


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


def generate_theoretical_ion_space(ms1_results_file, sites_file):
    return theoretical_glycopeptide.main(ms1_results_file, sites_file, None)


def match_deconvoluted_ions(theoretical_ion_space, deconvoluted_spectra, ms1_match_tolerance, ms2_match_tolerance):
    return match_ions.match_frags(theoretical_ion_space, deconvoluted_spectra, ms1_match_tolerance, ms2_match_tolerance, None)


def postprocess_matches(matched_ions_file):
    return postprocess.main(matched_ions_file)


def prepare_model_file(postprocessed_ions_file, method="default", rscript_path="Rscript", out=None):
    # FILE - Location of library, for non-bundled installed execution
    # STANDARD - The gold standard file to build a model froms
    # TEST - The file to test and classify
    # OUT - Final path to the result file
    # RSCRIPT - Path to the Rscript executable
    cmd_args = dict(FILE=windows_to_unix_path(os.path.dirname(os.path.abspath(__file__))),
                    MODEL=windows_to_unix_path(postprocessed_ions_file),
                    METHOD=method,
                    OUT="",
                    RSCRIPT=rscript_path)
    if(out is not None):
        cmd_args["OUT"] = "-o %s" % out
    cmd_str = """{RSCRIPT} --vanilla "{FILE}"/R/prepare_model.R -t "{MODEL}" -m {METHOD} {OUT}""".format(**cmd_args)
    print(cmd_str)
    cmd = subprocess.Popen(cmd_str, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Script writes result file to STDOUT
    stdout, stderr = cmd.communicate()
    # for stream in [stdout, stderr]:
    #     print(stream[:100])
    if cmd.returncode != 0:
        raise RscriptException(stderr)
    else:
        return stdout.split("\n")[-1]


def classify_data_by_gold_standard(postprocessed_ions_file, gold_standard_file, method='default', rscript_path="Rscript", out=None):
    # FILE - Location of library, for non-bundled installed execution
    # STANDARD - The gold standard file to build a model from
    # TEST - The file to test and classify
    # OUT - Final path to the result file
    # RSCRIPT - Path to the Rscript executable
    cmd_args = dict(FILE=windows_to_unix_path(os.path.dirname(os.path.abspath(__file__))),
                    STANDARD=windows_to_unix_path(gold_standard_file),
                    TEST=windows_to_unix_path(postprocessed_ions_file),
                    OUT="",
                    METHOD=method,
                    RSCRIPT=rscript_path)
    if(out is not None):
        cmd_args["OUT"] = "-o %s" % out
    cmd_str = """{RSCRIPT} --vanilla "{FILE}"/R/rf_classify.R -g "{STANDARD}" -t "{TEST}" -m {METHOD} {OUT}""".format(**cmd_args)
    print(cmd_str)
    cmd = subprocess.Popen(cmd_str, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Script writes result file to STDOUT
    stdout, stderr = cmd.communicate()
    # for stream in [stdout, stderr]:
    #     print(stream[:100])
    if cmd.returncode != 0:
        raise RscriptException(stderr)
    else:
        return stdout.split("\n")[-1]


def windows_to_unix_path(path):
    return path.replace("\\", "/")


class RscriptException(Exception):
    pass


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


#atexit.register(lambda: clean_up_files(*intermediary_files))
# atexit.register(__intermediary_files)


def main(ms1_results_file, glycosylation_sites_file, deconvoluted_spectra_file,
         ms1_match_tolerance, ms2_match_tolerance, constant_modification_list=None,
         variable_modification_list=None,
         method="default", rscript_path="Rscript", gold_standard_file=None, out=None):
    theoretical_ion_space_file = generate_theoretical_ion_space(ms1_results_file, glycosylation_sites_file)
    # print(theoretical_ion_space_file)
    intermediary_files.append(theoretical_ion_space_file)

    matched_ions_file = match_deconvoluted_ions(theoretical_ion_space_file, deconvoluted_spectra_file, ms1_match_tolerance, ms2_match_tolerance)
    # print(matched_ions_file)
    intermediary_files.append(matched_ions_file)

    postprocessed_ions_file = postprocess_matches(matched_ions_file)
    # print(postprocessed_ions_file)
    intermediary_files.append(postprocessed_ions_file)

    # If there is no gold standard, then we generate a model that awaits labeling
    if(gold_standard_file is None):
        try:
            model_file = prepare_model_file(postprocessed_ions_file, method=method,
                                            rscript_path=rscript_path, out=out)
            print(model_file)
        except RscriptException, e:
            print(e)
            exit(-25)
    else:
        try:
            classification_results_file = classify_data_by_gold_standard(postprocessed_ions_file, method=method,
                                                                         gold_standard_file=gold_standard_file,
                                                                         rscript_path=rscript_path, out=out)
        except RscriptException, e:
            print(e)
            exit(-25)
        print(classification_results_file)
    exit(0)

if __name__ == '__main__':
    import argparse
    app = argparse.ArgumentParser()
    app.add_argument("--parameter-file", action="store", default=None)
    app.add_argument("-r", "--rscript-path", action="store", default="Rscript", help="Path to the Rscript executable to run the R component")
    app.add_argument("--ms1-results-file", action="store")
    app.add_argument("--glycosylation-sites-file", action="store")
    app.add_argument("--deconvoluted-spectra-file", action="store")
    app.add_argument("--method", action="store", default="default", help="Select the model method to use for classification")
    app.add_argument("--gold-standard-file", action="store", default=None, required=False)
    app.add_argument("--ms1-match-tolerance", type=float, action="store", default=match_ions.ms1_tolerance_default, help="Mass Error Tolerance for matching MS1 masses in PPM")
    app.add_argument("--ms2-match-tolerance", type=float, action="store", default=match_ions.ms2_tolerance_default, help="Mass Error Tolerance for matching MS2 masses in PPM")
    app.add_argument("--constant-modification-list", type=str, action="append", default=None, help="Pass the list of constant modifications to include in the sequence search space")
    app.add_argument("--variable-modification-list", type=str, action="append", default=None, help="Pass the list of variable modifications to include in the sequence search space")

    app.add_argument("--out", action="store", default=None)
    args = app.parse_args().__dict__
    if args['constant_modification_list'] is not None:
        args['constant_modification_list'] = map(uri_decode, args['constant_modification_list'])

    if args['variable_modification_list'] is not None:
        args['variable_modification_list'] = map(uri_decode, args['variable_modification_list'])

    param_file = args.pop("parameter_file", None)
    if param_file is not None:
        params = load_parameters(param_file)
        args.update(params)
    print(json.dumps(args, indent=4))
    main(**args)
