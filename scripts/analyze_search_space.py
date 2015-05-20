from glycresoft_ms2_classification.entry_point import (match_deconvoluted_ions, postprocess_matches,
                                                       classify_data_by_model, calculate_false_discovery_rate,
                                                       match_ions2, ClassifyTargetWithModelTask)

import argparse

app = argparse.ArgumentParser()

app.add_argument("--search-space", action="store", required=True)

app.add_argument(
    "--deconvoluted-spectra-file", action="store", required=True)
app.add_argument(
    "--method", action="store", default="full_random_forest",
    choices=set(ClassifyTargetWithModelTask.method_table),
    help="Select the model method to use for classification")
app.add_argument("--model-file", action="store", default="naive")
app.add_argument(
    "--ms1-match-tolerance", type=float, action="store",
    default=match_ions2.ms1_tolerance_default,
    help="Mass Error Tolerance for matching MS1 masses in PPM")
app.add_argument(
    "--ms2-match-tolerance", type=float, action="store",
    default=match_ions2.ms2_tolerance_default,
    help="Mass Error Tolerance for matching MS2 masses in PPM")
app.add_argument("--decoy-to-real-ratio", action="store", default=1, type=int, help="Number of\
    decoys per prediction sequence")
app.add_argument("--prefix-length", default=0, required=False, type=int,
                                     help="Length of peptide prefix to preserve when generating\
                                      random glycopeptides by shuffling.")
app.add_argument("--suffix-length", default=1, required=False, type=int,
                                     help="Length of peptide suffix to preserve when generating\
                                      random glycopeptides by shuffling.")


def main():
    args = app.parse_args()
    out = match_deconvoluted_ions(
        args.search_space,
        args.deconvoluted_spectra_file,
        args.ms1_match_tolerance,
        args.ms2_match_tolerance,
        4)
    out = postprocess_matches(out)
    out = classify_data_by_model(out, args.model_file)
    out = calculate_false_discovery_rate(out, deconvoluted_spectra=args.deconvoluted_spectra_file, model_file_path=args.model_file)



if __name__ == '__main__':
    main()
