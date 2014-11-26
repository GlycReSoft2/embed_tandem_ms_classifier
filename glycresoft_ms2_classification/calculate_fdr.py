import csv
import itertools
import random
import re
import functools

import classify_matches
import match_ions
import postprocess

from .structure.sequence import Sequence
from .structure.sequence import strip_modifications
from .structure.sequence import sequence_tokenizer_respect_sequons
from .structure.sequence import list_to_sequence
from .structure.stub_glycopeptides import StubGlycopeptide
from .structure import modification


def make_predicate(**kwargs):
    fn = functools.partial(predicate_base, **kwargs)
    setattr(fn, "sig", str(kwargs))
    return fn


def and_predicate(predicate, other):
    fn = lambda x: predicate(x) and other(x)
    setattr(fn, "sig", "({0}) and ({1})".format(predicate.sig, other.sig))
    return fn


def or_predicate(predicate, other):
    fn = lambda x: predicate(x) or other(x)
    setattr(fn, "sig", "({0}) and ({1})".format(predicate.sig, other.sig))
    return fn


def negate_predicate(predicate):
    fn = lambda x: not predicate(x)
    setattr(fn, "sig", "not ({0})".format(predicate.sig))
    return fn


def predicate_base(x, MS2_Score=0, meanCoverage=0, meanHexNAcCoverage=0, percentUncovered=1, MS1_Score=0,
                   peptideLens=0, Obs_Mass=0, numStubs=-1):
    return (x.MS2_Score >= MS2_Score) & (x.meanCoverage >= meanCoverage) & (x.percentUncovered <= percentUncovered) &\
           (x.MS1_Score >= MS1_Score) & (x.peptideLens >= peptideLens) & (x.Obs_Mass >= Obs_Mass) &\
           (x.numStubs >= numStubs) & (x.meanHexNAcCoverage >= meanHexNAcCoverage)


def build_shuffle_sequences(scored_matches, count=20, prefix_len=0, suffix_len=0, iter_max=None):
    '''
        :type scored_matches: pd.DataFrame
        :param scored_matches: pd.DataFrame having scored MS2 matches
        :param count: int > 0, number of random peptides to generate per match
        :param prefix_len: int >= 0, number of AA to preserve order on of the
                           start of the match sequence
        :param suffix_len: int >= 0, number of AA to preserve order on of the
                           end of the match sequence
    '''
    iter_max = count * 100 if iter_max is None else iter_max
    shuffles = list()
    for ind, row in scored_matches.iterrows():
        solutions = set()
        seq = Sequence(row.Glycopeptide_identifier).get_sequence()
        solutions.add(seq)
        iter_count = 0
        clone = sequence_tokenizer_respect_sequons(row.Glycopeptide_identifier)
        pref = clone[:prefix_len]
        suf = clone[-suffix_len:]
        body = clone[prefix_len:len(clone) - suffix_len]
        while(len(solutions) - 1 < count and iter_count < iter_max):
            for shuffle in itertools.permutations(body):
                clone = pref + list(shuffle) + suf
                if random.random() < (1.0/(len(solutions) + 1.0)):
                    solutions.add(list_to_sequence(clone))
                if(len(solutions) - 1) >= count:
                    break
            iter_count += 1
        solutions.discard(seq)
        decoys = []
        for shuffle in solutions:
            d = row.copy()
            d.Glycopeptide_identifier = shuffle + row.Glycan
            d._batch_id = ind
            decoys.append(d)
        shuffles.append(decoys)
    return shuffles


def random_glycopeptide_to_fragments(sequence_record):
    seq_obj = Sequence(sequence_record.Glycopeptide_identifier)
    glycan_map = {}
    modifications = []
    for i, (aa, mods) in enumerate(seq_obj):
        for mod in mods:
            if mod.name in {"Glycan", "HexNAc"}:
                glycan_map[i] = mod.name
            else:
                # Construct the set of acceptable reasons why this modification is here.
                # Infer the least permissive modification rule.
                why = mod.why_valid(aa, i)
                modifications.append(modification.Modification(why, (i,)))

    # Remove glycans from the sequence string to conform to the SequenceSpace expectations
    for site, glycan in glycan_map.items():
        # Don't discard anonymous HexNAcs. Downstream functions can handle them
        if glycan != "HexNAc":
            seq_obj.drop_modification(site, glycan)
    glycosylation_sites = glycan_map.keys()
    if sequence_record.Glycan is None:
        # Build the semicolon separated string for glycan compositions
        glycan_composition = []
        glycan_composition = [map(int, glycan.replace("Glycan", '').replace("[", "").replace("]", "").split(";"))
                              for glycan in glycan_map.values()]
        glycan_composition = map(sum, zip(*glycan_composition))
        glycan_composition_restring = "[" + ";".join(map(str, glycan_composition)) + "]"
    else:
        glycan_composition_restring = sequence_record.Glycan
    # Begin generating fragment ions
    b_type = seq_obj.get_fragments('B')
    b_ions = []
    b_ions_HexNAc = []
    for b in b_type:
        for fm in b:
            key = fm.get_fragment_name()
            if key == "B1" or re.search(r'B1\+', key):
                # B1 Ions aren't actually seen in reality, but are an artefact of the generation process
                # so do not include them in the output
                continue
            mass = fm.get_mass()
            if "HexNAc" in key:
                b_ions_HexNAc.append({"key": key, "mass": mass})
            else:
                b_ions.append({"key": key, "mass": mass})

    y_type = seq_obj.get_fragments('Y')
    y_ions = []
    y_ions_HexNAc = []
    for y in y_type:
        for fm in y:
            key = fm.get_fragment_name()
            mass = fm.get_mass()
            if "HexNAc" in key:
                y_ions_HexNAc.append({"key": key, "mass": mass})
            else:
                y_ions.append({"key": key, "mass": mass})

    peptide_seq = strip_modifications(seq_obj.get_sequence())
    pep_stubs = StubGlycopeptide(peptide_seq, None,
                                 len(glycosylation_sites), glycan_composition_restring)
    stub_ions = pep_stubs.get_stubs()
    oxonium_ions = pep_stubs.get_oxonium_ions()
    ions = {
        "MS1_Score": sequence_record.MS1_Score,
        "Obs_Mass": sequence_record.Obs_Mass,
        "Calc_mass": sequence_record.Calc_mass,
        "ppm_error": sequence_record.ppm_error,
        "Peptide": peptide_seq,
        "Peptide_mod": sequence_record.Peptide_mod,
        "Glycan": glycan_composition_restring,
        "vol": sequence_record.vol,
        "glyco_sites": len(glycan_map),
        "startAA": None,
        "endAA": None,
        "Seq_with_mod": seq_obj.get_sequence(),
        "bare_b_ions": b_ions,
        "b_ions_with_HexNAc": b_ions_HexNAc,
        "bare_y_ions": y_ions,
        "y_ions_with_HexNAc": y_ions_HexNAc,
        "pep_stub_ions": stub_ions,
        "Oxonium_ions": oxonium_ions,
        "Glycopeptide_identifier": seq_obj.get_sequence() + glycan_composition_restring,
        "_batch_id": sequence_record._batch_id
    }
    return ions


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
        "predicted_database": None,
        "decoy_database": None,
        "searched_spectra": None,
        "threshold": (getattr(predicate_fn, "sig", None))
    }
    return results_obj


def generate_decoy_match_results(scored_matches_path, decon_data, model_file_path,
                                 prefix_len=0, suffix_len=0, ms1_tolerance=1e-5,
                                 ms2_tolerance=2e-5, num_decoys_per_real_mass=20.0,
                                 method="full_random_forest",
                                 method_init_args=None, method_fit_args=None):
    scored_matches_frame = classify_matches.prepare_model_file(
        scored_matches_path)
    decoy_file_name = scored_matches_path[:-4] + ".decoy.ion_space.csv"

    random_sequences = build_shuffle_sequences(scored_matches_frame, int(num_decoys_per_real_mass),
                                               prefix_len=prefix_len, suffix_len=suffix_len)
    fh = open(decoy_file_name, "wb")

    fragments = (random_glycopeptide_to_fragments(decoy) for decoy in itertools.chain.from_iterable(random_sequences))

    # Fetch the first fragment set for headers
    frag_sample = fragments.next()
    writer = csv.DictWriter(fh, frag_sample.keys())
    # Write out fragments
    writer.writeheader()
    writer.writerow(frag_sample)
    writer.writerows(fragments)
    fh.close()
    print(decoy_file_name)
    match_file = match_ions.match_frags(
        decoy_file_name, decon_data, ms1_tolerance, ms2_tolerance)
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
         predicate_fns=(make_predicate(MS2_Score=0.5),), prefix_len=0, suffix_len=0, by_mod_sig=False,
         ms1_tolerance=1e-5, ms2_tolerance=2e-5, method="full_random_forest", method_init_args=None,
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
    scored_matches_frame = classify_matches.prepare_model_file(
        scored_matches_path)
    decoy_matches_frame = classify_matches.prepare_model_file(
        decoy_matches_path)
    results = []
    for predicate in predicate_fns:
        # Modification Group Specific
        if by_mod_sig:
            for mod_signature in set(scored_matches_frame.modificationSignature):
                scored = scored_matches_frame.modificationSignature == mod_signature
                decoy = decoy_matches_frame.modificationSignature == mod_signature
                results_obj = calculate_fdr(
                    scored_matches_frame.ix[scored], decoy_matches_frame.ix[decoy],
                    predicate, num_decoys_per_real_mass)
                if results_obj is None:
                    continue
                results_obj['predicted_database'] = scored_matches_path
                results_obj['decoy_database'] = decoy_matches_path
                results_obj['modificationSignature'] = mod_signature
                results.append(results_obj)
        # Overall
        results_obj = calculate_fdr(
            scored_matches_frame, decoy_matches_frame,
            predicate, num_decoys_per_real_mass)
        if results_obj is None:
            continue
        results_obj['predicted_database'] = scored_matches_path
        results_obj['decoy_database'] = decoy_matches_path
        results_obj['modificationSignature'] = "*"
        results.append(results_obj)

    results_fh = open(outfile_path, 'wb')
    results_writer = csv.DictWriter(results_fh, results[0].keys())
    results_writer.writeheader()
    results_writer.writerows(results)
    results_fh.close()
    return results


def default_predicates():
    predicates = [make_predicate(MS2_Score=i)
                  for i in [0.2, 0.4, 0.6, 0.8]]
    predicates.extend([make_predicate(MS2_Score=i, numStubs=1)
                       for i in [0.2, 0.4, 0.6, 0.8]])
    predicates.extend([make_predicate(MS2_Score=i, peptideLens=10)
                       for i in [0.2, 0.4, 0.6, 0.8]])
    predicates.extend([and_predicate(make_predicate(MS2_Score=i),
                       negate_predicate(make_predicate(peptideLens=10)))
                       for i in [0.2, 0.4, 0.6, 0.8]])
    predicates.extend([make_predicate(
        MS2_Score=i, peptideLens=10, numStubs=1) for i in [0.2, 0.4, 0.6, 0.8]])
    predicates.extend([make_predicate(
        MS2_Score=i, numStubs=1, meanHexNAcCoverage=.5) for i in [0.2, 0.4, 0.6, 0.8]])

    return predicates


if __name__ == '__main__':
    import argparse
    app = argparse.ArgumentParser()
    app.add_argument(
        "--decon_data", default=None, help="MS2 deconvolution output file path")
    app.add_argument("--model_file", default=None,
                     help="The model to train a classifier on for scoring decoys")
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
