import os
import sys
import subprocess

# Ensure our scripts are imported before any others since we can't depend upon
# the packaging middleware to handle relative imports, or should we vendorize
# pip or easy_install?
sys.path = [(os.path.dirname(os.path.abspath(__file__)))] + sys.path

# Actions to be bound to model
import theoretical_glycopeptide
import match_ions
import postprocess


def generate_theoretical_ion_space(ms1_results_file, sites_file):
    return theoretical_glycopeptide.main(ms1_results_file, sites_file, None)


def match_deconvoluted_ions(theoretical_ion_space, deconvoluted_spectra):
    return match_ions.match_frags(theoretical_ion_space, deconvoluted_spectra, None)


def postprocess_matches(matched_ions_file):
    return postprocess.main(matched_ions_file)


def classify_data_by_gold_standard(postprocessed_ions_file, gold_standard_file, rscript_path="Rscript", out=None):
    # FILE - Location of library, for non-bundled installed execution
    # STANDARD - The gold standard file to build a model from
    # TEST - The file to test and classify
    # OUT - Final path to the result file
    # RSCRIPT - Path to the Rscript executable
    cmd_args = dict(FILE=windows_to_unix_path(os.path.dirname(os.path.abspath(__file__))),
                    STANDARD=windows_to_unix_path(gold_standard_file),
                    TEST=windows_to_unix_path(postprocessed_ions_file),
                    OUT="",
                    RSCRIPT=rscript_path)
    if(out is not None):
        cmd_args["OUT"] = "-o %s" % out
    cmd_str = """{RSCRIPT} --vanilla "{FILE}"/R/rf_classify.R -g "{STANDARD}" -t "{TEST}" {OUT}""".format(**cmd_args)
    cmd = subprocess.Popen(cmd_str,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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


def main(ms1_results_file, glycosylation_sites_file, deconvoluted_spectra_file, rscript_path="Rscript", gold_standard_file=None, out=None):
    theoretical_ion_space_file = generate_theoretical_ion_space(ms1_results_file, glycosylation_sites_file)
    print(theoretical_ion_space_file)
    matched_ions_file = match_deconvoluted_ions(theoretical_ion_space_file, deconvoluted_spectra_file)
    print(matched_ions_file)
    postprocessed_ions_file = postprocess_matches(matched_ions_file)
    print(postprocessed_ions_file)
    try:
        classification_results_file = classify_data_by_gold_standard(postprocessed_ions_file, gold_standard_file, rscript_path, out)
    except RscriptException, e:
        print(e)
        exit(-25)
    print(classification_results_file)
    exit(0)

if __name__ == '__main__':
    import argparse
    app = argparse.ArgumentParser()
    app.add_argument("-r", "--rscript-path", action="store", default="Rscript", help="Path to the Rscript executable to run the R component")
    app.add_argument("--ms1-results-file", action="store")
    app.add_argument("--glycosylation-sites-file", action="store")
    app.add_argument("--deconvoluted-spectra-file", action="store")
    app.add_argument("--gold-standard-file", action="store", default=None)
    app.add_argument("--out", action="store", default=None)
    args = app.parse_args()
    main(**args.__dict__)
