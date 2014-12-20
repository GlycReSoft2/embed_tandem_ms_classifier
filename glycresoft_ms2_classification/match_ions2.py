import json
import math
from os.path import splitext
from collections import Counter, defaultdict

import multiprocessing
import functools
import itertools

import logging
multiprocessing.log_to_stderr()

from glycresoft_ms2_classification.error_code_interface import NoIonsMatchedException
from glycresoft_ms2_classification.utils import try_deserialize, try_get_outfile, collectiontools
from glycresoft_ms2_classification.data_io.bupid_topdown_deconvoluter import BUPIDYamlParser
from glycresoft_ms2_classification.data_io import default_loader, default_serializer
from glycresoft_ms2_classification.data_io import DeconIOBase
from glycresoft_ms2_classification.data_io import ParallelParser

Proton = 1.007276035

ms1_tolerance_default = 10e-6
ms2_tolerance_default = 20e-6

MACHINE_ACQUISITION_RANGE = 3200 * 24


def save_interm(data, name):
    json.dump(data, open(name + ".json", 'wb'))


def split_decon_data_by_index(decon_data, splitting_index):
    '''
    It might be useful to know which scans in the data matched and which did not.
    :param decon_data: dict-like yaml data describing the deconvoluted spectra we matched against
    :param splitting_index: dict-like record of which peaks matched
    '''
    matched = {}
    no_matched = {}
    for i, scan in decon_data:
        count = splitting_index[i]
        if count != 0:
            matched[i] = scan
        else:
            no_matched[i] = scan

    logging.info("Matched {matched} spectra".format(matched=len(matched)))
    logging.info(
        "Did not match {unmatched} spectra".format(unmatched=len(no_matched)))
    wrapped_matched = DeconIOBase.prepare(matched)
    wrapped_no_matched = DeconIOBase.prepare(no_matched)
    return wrapped_matched, wrapped_no_matched


def save_split_decon_data(matched, unmatched, file_stem):
    default_serializer(matched, file_stem + "-matched.pkl")
    default_serializer(unmatched, file_stem + "-unmatched.pkl")


decon_format_lookup = {
    "bupid_yaml": BUPIDYamlParser,
    "default": default_loader
}


def match_observed_to_theoretical(theoretical, observed_ions, ms1_tolerance, ms2_tolerance):
    obs_mass = float(theoretical['Calc_mass'])
    results = []
    did_match_counter = Counter()
    fabs = math.fabs
    for tandem_ms_ind, real_ions in observed_ions:
        scan_id = real_ions.scan_ids[0]

        precursor_ppm = (
            (real_ions.neutral_mass - obs_mass) / real_ions.neutral_mass)

        if fabs(precursor_ppm) <= ms1_tolerance:
            did_match_counter[tandem_ms_ind] += 1
            Oxonium_ions = theoretical['Oxonium_ions']

            oxoniums = []
            for ox_ion in Oxonium_ions:
                for ob_ions in real_ions.tandem_data:
                    oxonium_ppm = (
                        ((ob_ions.mass + Proton) - ox_ion) / (ob_ions.mass + Proton))
                    if fabs(oxonium_ppm) <= ms2_tolerance:
                        oxoniums.append(
                            {"ion": (ob_ions.mass + Proton), "ppm_error": oxonium_ppm * 1e6,
                             "key": ox_ion, "scan_id": tandem_ms_ind})

            # checking for b and y ions and ions with HexNAc:
            b_ions = theoretical['bare_b_ions']
            b_len = float(len(b_ions))
            b_type = []
            all_b_ions = []
            for theo_ions in b_ions:
                frag_mass = float(theo_ions["mass"])
                prot_ion = frag_mass - Proton
                for obs_ions in real_ions.tandem_data:
                    tandem_ppm = float(
                        (obs_ions.mass - (prot_ion)) / obs_ions.mass)
                    if fabs(tandem_ppm) <= ms2_tolerance:
                        b_type.append(
                            {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
                        all_b_ions.append(
                            {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})

            b_HexNAc_ions = theoretical['b_ions_with_HexNAc']
            b_HexNAc_type = []
            for theo_ions in b_HexNAc_ions:
                frag_mass = float(theo_ions["mass"])
                prot_ion = frag_mass - Proton
                for obs_ions in real_ions.tandem_data:
                    tandem_ppm = float(
                        (obs_ions.mass - (prot_ion)) / obs_ions.mass)
                    if fabs(tandem_ppm) <= ms2_tolerance:
                        b_HexNAc_type.append(
                            {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
                        all_b_ions.append(
                            {"obs_ion": (obs_ions.mass + Proton),
                             "key": theo_ions["key"].split("+")[0],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})

            y_ions = list(theoretical['bare_y_ions'])
            y_len = float(len(y_ions))
            y_type = []
            all_y_ions = []
            for theo_ions in y_ions:
                frag_mass = float(theo_ions["mass"])
                prot_ion = frag_mass - Proton
                for obs_ions in real_ions.tandem_data:
                    tandem_ppm = float(
                        (obs_ions.mass - (prot_ion)) / obs_ions.mass)
                    if fabs(tandem_ppm) <= ms2_tolerance:
                        y_type.append(
                            {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
                        all_y_ions.append(
                            {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})

            y_HexNAc_ions = theoretical['y_ions_with_HexNAc']
            y_HexNAc_type = []
            for theo_ions in y_HexNAc_ions:
                frag_mass = float(theo_ions["mass"])
                prot_ion = frag_mass - Proton
                for obs_ions in real_ions.tandem_data:
                    tandem_ppm = float(
                        (obs_ions.mass - (prot_ion)) / obs_ions.mass)
                    if fabs(tandem_ppm) <= ms2_tolerance:
                        y_HexNAc_type.append(
                            {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
                        all_y_ions.append(
                            {"obs_ion": (obs_ions.mass + Proton),
                             "key": theo_ions["key"].split("+")[0],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})

            # checking for stub ions
            stub_ions = theoretical['pep_stub_ions']
            stub_type = []
            for theo_ions in stub_ions:
                frag_mass = float(theo_ions["mass"])
                prot_ion = frag_mass - Proton
                for obs_ions in real_ions.tandem_data:
                    tandem_ppm = float(
                        (obs_ions.mass - (prot_ion)) / obs_ions.mass)
                    if fabs(tandem_ppm) <= ms2_tolerance:
                        stub_type.append(
                            {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
            results.append({
                "MS1_Score": theoretical["MS1_Score"], "Obs_Mass": theoretical["Obs_Mass"],
                "Calc_mass": theoretical["Calc_mass"],
                "ppm_error": theoretical["ppm_error"], "Peptide": theoretical["Peptide"],
                "Peptide_mod": theoretical["Peptide_mod"],
                "Glycan": theoretical["Glycan"],
                "vol": theoretical["vol"], "glyco_sites": theoretical["glyco_sites"],
                "startAA": theoretical["startAA"], "endAA": theoretical["endAA"],
                "Seq_with_mod": theoretical["Seq_with_mod"],
                "Glycopeptide_identifier": theoretical["Seq_with_mod"] + theoretical["Glycan"],
                "Oxonium_ions": oxoniums, "bare_b_ions": b_type, "possible_b_ions_HexNAc": len(b_HexNAc_ions),
                "total_b_ions": b_len, "bare_y_ions": y_type, "possible_y_ions_HexNAc": len(y_HexNAc_ions),
                "total_y_ions": y_len, "b_ions_with_HexNAc": b_HexNAc_type,
                "y_ions_with_HexNAc": y_HexNAc_type,
                "b_ion_coverage": all_b_ions, "y_ion_coverage": all_y_ions, "Stub_ions": stub_type,
                "scan_id": scan_id,
                "scan_id_range": []
            })

    return results, did_match_counter


def merge_by_glycopeptide_sequence(matches):
    groups = collectiontools.groupby(matches, lambda x: x["Glycopeptide_identifier"])
    merged_matches = []
    for glycopeptide, matches in groups.items():
        merged_matches.append(merge_matches(matches))

    return merged_matches

fields_to_merge = [
    "Oxonium_ions",
    "bare_b_ions",
    "bare_y_ions",
    "b_ions_with_HexNAc",
    "y_ions_with_HexNAc",
    "b_ion_coverage",
    "y_ion_coverage",
    "Stub_ions"
]


def merge_matches(matches):
    merged = {}
    matches = list(matches)
    scan_ids = set()
    for match in matches:
        scan_ids.add(match["scan_id"])
    for k, v in matches[0].items():
        if k in fields_to_merge:
            merged[k] = merge_ion_matches([m[k] for m in matches])
        else:
            merged[k] = v
    merged["scan_id_range"] = list(scan_ids)
    return merged


def merge_ion_matches(matches):
    groups = collectiontools.groupby(itertools.chain.from_iterable(matches),
                                     lambda x: x["key"])
    best_matches = []
    fabs = math.fabs
    # for match in itertools.chain.from_iterable(matches):
    #     groups[match["key"]].append(match)
    for key, matched_key in groups.items():
        best_match = matched_key[0]
        best_ppm = fabs(best_match["ppm_error"])
        ppm_to_scan_id = {best_match["scan_id"]: best_match["ppm_error"]}
        for match in matched_key[1:]:
            ppm_to_scan_id[match["scan_id"]] = match["ppm_error"]
            if fabs(match["ppm_error"]) < best_ppm:
                best_match = match
                best_ppm = fabs(match["ppm_error"])
        best_match["scan_map"] = ppm_to_scan_id
        best_matches.append(best_match)
    return best_matches


def match_frags(theo_fragment_file, decon_data, ms1_tolerance=ms1_tolerance_default,
                ms2_tolerance=ms2_tolerance_default, split_decon_data=False,
                tag="", decon_data_format='bupid_yaml', n_processes=4,
                outfile=None):
    '''
    :param theo_fragment_file: path to file containing all theoretical peptides
    :param decon_data: path to yaml file from earlier MS1 analysis
    :param ms1_tolerance: float tolerance for matching intact mass
    :param ms2_tolerance: float tolerance for matching fragment ions
    :param split_decon_data: bool flag whether to emit a file for all scans without matches
                             and a file for all scans with matches.
    :param tag: name to associated with match
    :param outfile: path to the file to write output to. Defaults to theo_fragment_file + ".match_frags".
    '''

    pool = multiprocessing.Pool(n_processes)

    data = ParallelParser(decon_format_lookup[decon_data_format], (decon_data,))

    theoretical_search_space = try_deserialize(theo_fragment_file)
    if outfile is None:
        outfile = try_get_outfile(theoretical_search_space, "match_frags")

    metadata = theoretical_search_space["metadata"]
    metadata["ms1_ppm_tolerance"] = ms1_tolerance
    metadata["ms2_ppm_tolerance"] = ms2_tolerance

    if "tag" not in metadata:
        metadata["tag"] = tag
    if metadata["tag"] == "decoy":
        tag = "decoy"
    elif tag is None:
        tag = metadata["tag"]

    theoretical_fragments = theoretical_search_space[
        "theoretical_search_space"]

    print(len(theoretical_fragments))

    results = []
    did_match_counter = Counter()

    data = data.await()

    task_fn = functools.partial(match_observed_to_theoretical, ms1_tolerance=ms1_tolerance,
                                ms2_tolerance=ms2_tolerance, observed_ions=data)
    matching_process = pool.imap_unordered(task_fn, theoretical_fragments,
                                           chunksize=len(theoretical_fragments)/n_processes)
    for matches, counter in matching_process:
        results.extend(matches)
        did_match_counter += counter

    pool.close()
    pool.join()
    if(len(results) < 1):
        raise NoIonsMatchedException(
            "No matches found from theoretical ions in MS2 deconvoluted results")

    # save_unmerged = multiprocessing.Process(target=save_interm, args=(results, "unmerged-matches"))
    # save_unmerged.start()

    merged_results = merge_by_glycopeptide_sequence(results)

    results = {
        "metadata": metadata,
        "matched_ions": merged_results
    }
    f = open(outfile + '.json', 'wb')
    json.dump(results, f)
    f.close()

    if(split_decon_data):
        match, no_match = split_decon_data_by_index(data, did_match_counter)
        save_split_decon_data(match, no_match,
                              "{root}.{tag}".format(root=splitext(decon_data)[0],
                                                    tag=tag))

    return outfile + '.json', results


if __name__ == '__main__':
    import argparse
    app = argparse.ArgumentParser("match-frags")
    app.add_argument("theoretical_fragments")
    app.add_argument("observed_ions")
    app.add_argument("-s", "--split-observations",
                     action="store_true", default=False, required=False)
    app.add_argument("-f", "--format", action="store", default="bupid_yaml",
                     choices=decon_format_lookup.keys(), required=False)
    app.add_argument("-t", "--tag", action="store", default="", required=False)

    args = app.parse_args()

    match_frags(args.theoretical_fragments, args.observed_ions, tag=args.tag,
                decon_data_format=args.format, split_decon_data=args.split_observations)
