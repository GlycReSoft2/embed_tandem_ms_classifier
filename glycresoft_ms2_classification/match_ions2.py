import json
import math
import multiprocessing
import functools
import itertools
import logging
import random
try:
    logger = logging.getLogger(__name__)
except:
    pass
from os.path import splitext
from collections import Counter, defaultdict

import sqlitedict

from glycresoft_ms2_classification.error_code_interface import NoIonsMatchedException
from glycresoft_ms2_classification.utils import collectiontools
from glycresoft_ms2_classification.ms.bupid_topdown_deconvoluter import BUPIDYamlParser
from glycresoft_ms2_classification.ms import default_loader, default_serializer
from glycresoft_ms2_classification.ms import DeconIOBase, MSMSSqlDB, ObservedPrecursorSpectrum
from glycresoft_ms2_classification.utils.parallel_opener import ParallelParser

use_cython = True
try:
    from glycresoft_ms2_classification.ms import ion_matching
    c_match_observed_to_theoretical_sql = ion_matching.match_observed_wrap
except:
    use_cython = False

Proton = 1.007276035

ms1_tolerance_default = 10e-6
ms2_tolerance_default = 20e-6


def save_interm(data, name):
    json.dump(data, open("./{0}.json".format(name), 'wb'))


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

    logger.info("Matched {matched} spectra".format(matched=len(matched)))
    logger.info(
        "Did not match {unmatched} spectra".format(unmatched=len(no_matched)))
    wrapped_matched = DeconIOBase.prepare(matched)
    wrapped_no_matched = DeconIOBase.prepare(no_matched)
    return wrapped_matched, wrapped_no_matched


def save_split_decon_data(matched, unmatched, file_stem):
    logger.info("Saving spectra to %s", file_stem)
    default_serializer(matched, file_stem + "matched.pkl")
    default_serializer(unmatched, file_stem + "unmatched.pkl")


def combine_annotations(accumulator, collected):
    for k, v in collected.items():
        accumulator[k].extend(v)


def annotate_ions(observed_ions, annotations):
    for iterkey, precursor in observed_ions:
        for tandem_ion in precursor.tandem_data:
            tandem_ion.annotation = annotations[tandem_ion.id]


decon_format_lookup = {
    "bupid_yaml": BUPIDYamlParser,
    "db": MSMSSqlDB,
    "default": default_loader
}


def dispatcher(index_range, theoretical_space_path, worker_fn):
    range_start, range_stop = index_range
    theoretical_search_space = sqlitedict.SqliteDict(theoretical_space_path, tablename="theoretical_search_space")
    did_match_counter = Counter()
    annotate_accumulator = defaultdict(list)
    fragment_match_store = sqlitedict.SqliteDict(theoretical_space_path, tablename="matched_ions", journal_mode="OFF")
    for i in xrange(range_start, range_stop):
        theoretical = theoretical_search_space.get(i, None)
        if theoretical is None:
            continue
        matches, counter, annotater = worker_fn(theoretical)
        for m in matches:
            fragment_match_store[i] = m
            i += 1
        did_match_counter += counter
        combine_annotations(annotate_accumulator, annotater)
    return did_match_counter, annotate_accumulator


def match_observed_to_theoretical_sql(theoretical, observed_ions_conn_string, ms1_tolerance, ms2_tolerance):
    db = MSMSSqlDB(observed_ions_conn_string)
    theoretical_sequence = theoretical["Seq_with_mod"] + theoretical["Glycan"]
    obs_mass = float(theoretical['Calc_mass'])
    results = []
    did_match_counter = Counter()
    annotate = defaultdict(list)
    fabs = math.fabs

    scan_id_range = []

    oxoniums = []
    b_type = []
    all_b_ions = []
    b_HexNAc_type = []
    y_type = []
    all_y_ions = []
    y_HexNAc_type = []
    stub_type = []

    for real_ions in db.ppm_match_tolerance_search(obs_mass, ms1_tolerance):
        real_ions = ObservedPrecursorSpectrum.from_sql(real_ions, db)
        scan_id = real_ions.scan_ids[0]
        scan_id_range.append(scan_id)
        tandem_ms_ind = real_ions.id
        did_match_counter[tandem_ms_ind] += 1

        Oxonium_ions = theoretical['Oxonium_ions']

        for ox_ion in Oxonium_ions:
            try:
                kind = ox_ion['key']
                mass = ox_ion['mass']
            except:
                kind = mass = ox_ion
            for obs_ions in real_ions.tandem_data:
                oxonium_ppm = (
                    ((obs_ions.mass + Proton) - mass) / (mass))
                if fabs(oxonium_ppm) <= ms2_tolerance:
                    oxoniums.append(
                        {"ion": (obs_ions.mass + Proton), "ppm_error": oxonium_ppm * 1e6,
                         "key": kind, "scan_id": obs_ions.id})
                    annotate[obs_ions.id].append((theoretical_sequence, kind))

        # checking for b and y ions and ions with HexNAc:
        b_ions = theoretical['bare_b_ions']
        b_len = float(len(b_ions))
        for theo_ions in b_ions:
            frag_mass = float(theo_ions["mass"])
            deprotonated_ion = frag_mass - Proton
            for obs_ions in real_ions.tandem_data:
                tandem_ppm = float(
                    (obs_ions.mass - (deprotonated_ion)) / deprotonated_ion)
                if fabs(tandem_ppm) <= ms2_tolerance:
                    b_type.append(
                        {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                         "ppm_error": tandem_ppm * 1e6, "scan_id": obs_ions.id})
                    all_b_ions.append(
                        {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                         "ppm_error": tandem_ppm * 1e6, "scan_id": obs_ions.id})
                    annotate[obs_ions.id].append((theoretical_sequence, theo_ions["key"]))

        b_HexNAc_ions = theoretical['b_ions_with_HexNAc']
        for theo_ions in b_HexNAc_ions:
            frag_mass = float(theo_ions["mass"])
            deprotonated_ion = frag_mass - Proton
            for obs_ions in real_ions.tandem_data:
                tandem_ppm = float(
                    (obs_ions.mass - (deprotonated_ion)) / deprotonated_ion)
                if fabs(tandem_ppm) <= ms2_tolerance:
                    b_HexNAc_type.append(
                        {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                         "ppm_error": tandem_ppm * 1e6, "scan_id": obs_ions.id})
                    all_b_ions.append(
                        {"obs_ion": (obs_ions.mass + Proton),
                         "key": theo_ions["key"].split("+")[0],
                         "ppm_error": tandem_ppm * 1e6, "scan_id": obs_ions.id})
                    annotate[obs_ions.id].append((theoretical_sequence, theo_ions["key"]))

        y_ions = list(theoretical['bare_y_ions'])
        y_len = float(len(y_ions))
        for theo_ions in y_ions:
            frag_mass = float(theo_ions["mass"])
            deprotonated_ion = frag_mass - Proton
            for obs_ions in real_ions.tandem_data:
                tandem_ppm = float(
                    (obs_ions.mass - (deprotonated_ion)) / deprotonated_ion)
                if fabs(tandem_ppm) <= ms2_tolerance:
                    y_type.append(
                        {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                         "ppm_error": tandem_ppm * 1e6, "scan_id": obs_ions.id})
                    all_y_ions.append(
                        {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                         "ppm_error": tandem_ppm * 1e6, "scan_id": obs_ions.id})
                    annotate[obs_ions.id].append((theoretical_sequence, theo_ions["key"]))

        y_HexNAc_ions = theoretical['y_ions_with_HexNAc']
        for theo_ions in y_HexNAc_ions:
            frag_mass = float(theo_ions["mass"])
            deprotonated_ion = frag_mass - Proton
            for obs_ions in real_ions.tandem_data:
                tandem_ppm = float(
                    (obs_ions.mass - (deprotonated_ion)) / deprotonated_ion)
                if fabs(tandem_ppm) <= ms2_tolerance:
                    y_HexNAc_type.append(
                        {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                         "ppm_error": tandem_ppm * 1e6, "scan_id": obs_ions.id})
                    all_y_ions.append(
                        {"obs_ion": (obs_ions.mass + Proton),
                         "key": theo_ions["key"].split("+")[0],
                         "ppm_error": tandem_ppm * 1e6, "scan_id": obs_ions.id})
                    annotate[obs_ions.id].append((theoretical_sequence, theo_ions["key"]))

        # checking for stub ions
        stub_ions = theoretical['pep_stub_ions']

        for theo_ions in stub_ions:
            frag_mass = float(theo_ions["mass"])
            deprotonated_ion = frag_mass - Proton
            for obs_ions in real_ions.tandem_data:
                tandem_ppm = float(
                    (obs_ions.mass - (deprotonated_ion)) / deprotonated_ion)
                if fabs(tandem_ppm) <= ms2_tolerance:
                    stub_type.append(
                        {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                         "ppm_error": tandem_ppm * 1e6, "scan_id": obs_ions.id})
                    annotate[obs_ions.id].append((theoretical_sequence, theo_ions["key"]))

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
        "possible_b_ions_HexNAc": len(b_HexNAc_ions),
        "total_b_ions": b_len, "possible_y_ions_HexNAc": len(y_HexNAc_ions),
        "total_y_ions": y_len,

        "Oxonium_ions": merge_ion_matches(oxoniums),
        "bare_b_ions": merge_ion_matches(b_type),
        "bare_y_ions": merge_ion_matches(y_type),
        "b_ions_with_HexNAc": merge_ion_matches(b_HexNAc_type),
        "y_ions_with_HexNAc": merge_ion_matches(y_HexNAc_type),
        "b_ion_coverage": merge_ion_matches(all_b_ions),
        "y_ion_coverage": merge_ion_matches(all_y_ions),
        "Stub_ions": merge_ion_matches(stub_type),
        "scan_id": scan_id_range[0],
        "scan_id_range": scan_id_range
    })

    return results, did_match_counter, annotate


def iterative_pooling_matcher(theoretical, observed_ions, ms1_tolerance, ms2_tolerance):
    theoretical_sequence = theoretical["Seq_with_mod"] + theoretical["Glycan"]
    calc_mass = float(theoretical['Calc_mass'])
    results = []
    did_match_counter = Counter()
    annotate = defaultdict(list)
    fabs = math.fabs
    for tandem_ms_ind, real_ions in observed_ions:
        scan_id = real_ions.scan_ids[0]

        precursor_ppm = (
            (real_ions.neutral_mass - calc_mass) / calc_mass)

        if fabs(precursor_ppm) <= ms1_tolerance:
            did_match_counter[tandem_ms_ind] += 1
            Oxonium_ions = theoretical['Oxonium_ions']

            oxoniums = []
            for ox_ion in Oxonium_ions:
                for obs_ions in real_ions.tandem_data:
                    oxonium_ppm = (
                        ((obs_ions.mass + Proton) - ox_ion) / ox_ion)
                    if fabs(oxonium_ppm) <= ms2_tolerance:
                        oxoniums.append(
                            {"ion": (obs_ions.mass + Proton), "ppm_error": oxonium_ppm * 1e6,
                             "key": ox_ion, "scan_id": tandem_ms_ind})
                        annotate[obs_ions.id].append((theoretical_sequence, ox_ion))

            # checking for b and y ions and ions with HexNAc:
            b_ions = theoretical['bare_b_ions']
            b_len = float(len(b_ions))
            b_type = []
            all_b_ions = []
            for theo_ions in b_ions:
                frag_mass = float(theo_ions["mass"])
                deprotonated_ion = frag_mass - Proton
                for obs_ions in real_ions.tandem_data:
                    tandem_ppm = float(
                        (obs_ions.mass - (deprotonated_ion)) / deprotonated_ion)
                    if fabs(tandem_ppm) <= ms2_tolerance:
                        b_type.append(
                            {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
                        all_b_ions.append(
                            {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
                        annotate[obs_ions.id].append((theoretical_sequence, theo_ions["key"]))

            b_HexNAc_ions = theoretical['b_ions_with_HexNAc']
            b_HexNAc_type = []
            for theo_ions in b_HexNAc_ions:
                frag_mass = float(theo_ions["mass"])
                deprotonated_ion = frag_mass - Proton
                for obs_ions in real_ions.tandem_data:
                    tandem_ppm = float(
                        (obs_ions.mass - (deprotonated_ion)) / deprotonated_ion)
                    if fabs(tandem_ppm) <= ms2_tolerance:
                        b_HexNAc_type.append(
                            {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
                        all_b_ions.append(
                            {"obs_ion": (obs_ions.mass + Proton),
                             "key": theo_ions["key"].split("+")[0],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
                        annotate[obs_ions.id].append((theoretical_sequence, theo_ions["key"]))

            y_ions = list(theoretical['bare_y_ions'])
            y_len = float(len(y_ions))
            y_type = []
            all_y_ions = []
            for theo_ions in y_ions:
                frag_mass = float(theo_ions["mass"])
                deprotonated_ion = frag_mass - Proton
                for obs_ions in real_ions.tandem_data:
                    tandem_ppm = float(
                        (obs_ions.mass - (deprotonated_ion)) / deprotonated_ion)
                    if fabs(tandem_ppm) <= ms2_tolerance:
                        y_type.append(
                            {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
                        all_y_ions.append(
                            {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
                        annotate[obs_ions.id].append((theoretical_sequence, theo_ions["key"]))

            y_HexNAc_ions = theoretical['y_ions_with_HexNAc']
            y_HexNAc_type = []
            for theo_ions in y_HexNAc_ions:
                frag_mass = float(theo_ions["mass"])
                deprotonated_ion = frag_mass - Proton
                for obs_ions in real_ions.tandem_data:
                    tandem_ppm = float(
                        (obs_ions.mass - (deprotonated_ion)) / deprotonated_ion)
                    if fabs(tandem_ppm) <= ms2_tolerance:
                        y_HexNAc_type.append(
                            {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
                        all_y_ions.append(
                            {"obs_ion": (obs_ions.mass + Proton),
                             "key": theo_ions["key"].split("+")[0],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
                        annotate[obs_ions.id].append((theoretical_sequence, theo_ions["key"]))

            # checking for stub ions
            stub_ions = theoretical['pep_stub_ions']
            stub_type = []
            for theo_ions in stub_ions:
                frag_mass = float(theo_ions["mass"])
                deprotonated_ion = frag_mass - Proton
                for obs_ions in real_ions.tandem_data:
                    tandem_ppm = float(
                        (obs_ions.mass - (deprotonated_ion)) / deprotonated_ion)
                    if fabs(tandem_ppm) <= ms2_tolerance:
                        stub_type.append(
                            {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                             "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
                        annotate[obs_ions.id].append((theoretical_sequence, theo_ions["key"]))

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

    return results, did_match_counter, annotate


def merge_by_glycopeptide_sequence(matches):
    groups = collectiontools.groupby(matches.values(), lambda x: x["Glycopeptide_identifier"], kind="sqlist")
    merged_matches = sqlitedict.SqliteDict(matches.filename, tablename="matched_ions")
    cntr = 0
    for glycopeptide, matches in groups.items():
        merged_matches[cntr] = merge_matches(matches)
        cntr += 1
    merged_matches.commit()
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
    for ion_type, ion_series in matches[0].items():
        if ion_type in fields_to_merge:
            merged[ion_type] = merge_ion_matches([m[ion_type] for m in matches])
        else:
            merged[ion_type] = ion_series
    merged["scan_id_range"] = list(scan_ids)
    return merged


def merge_ion_matches(matches):
    groups = collectiontools.groupby(itertools.chain.from_iterable(matches),
                                     lambda x: x["key"])
    best_matches = []
    fabs = math.fabs
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


def match_frags(db_file, decon_data, ms1_tolerance=ms1_tolerance_default,
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

    pool = multiprocessing.Pool(n_processes, maxtasksperchild=10)

    # theoretical_search_space = try_deserialize(theo_fragment_file)
    theoretical_search_space = sqlitedict.SqliteDict(
        db_file, tablename="theoretical_search_space")
    if isinstance(decon_data, str):
        if decon_data_format != "db" and splitext(decon_data)[1] != '.db':
            data = ParallelParser(decon_format_lookup[decon_data_format], (decon_data,))
            db = None
        else:
            data = None
            db = MSMSSqlDB(decon_data)
    elif isinstance(decon_data, MSMSSqlDB):
        data = None
        db = decon_data

    metadata = sqlitedict.SqliteDict(db_file, tablename="metadata")
    metadata["ms1_ppm_tolerance"] = ms1_tolerance
    metadata["ms2_ppm_tolerance"] = ms2_tolerance
    metadata["deconvoluted_spectra_file"] = decon_data
    if tag == "":
        tag = hex(random.randint(100000, 1000000))
    if "tag" not in metadata:
        metadata["tag"] = tag
    else:
        tag = metadata["tag"]
    metadata.commit()

    logging.info("Run tag: %s", tag)
    theoretical_fragments = theoretical_search_space.itervalues()
    fragment_match_store = sqlitedict.SqliteDict(db_file, tablename="matched_ions", journal_mode="OFF")

    did_match_counter = Counter()
    annotate_accumulator = defaultdict(list)

    if db is None and data is not None:
        data = data.await()
        logger.info("Indexing observed ions")
        db = data.to_db()
        logger.info("Index built.")
    elif db is not None:
        data = []

    match_fn = match_observed_to_theoretical_sql
    if use_cython:
        logger.info("Using Cython")
        match_fn = c_match_observed_to_theoretical_sql
    task_fn = functools.partial(match_fn, ms1_tolerance=ms1_tolerance,
                                ms2_tolerance=ms2_tolerance,
                                observed_ions_conn_string=db.connection_string)
    # Index counter for sqlitedict
    cntr = 0
    try:
        if n_processes > 1:
            logger.debug("Matching concurrently")
            matching_process = pool.imap_unordered(task_fn, theoretical_fragments,
                                                   chunksize=750)
            for matches, counter, annotater in matching_process:
                for m in matches:
                    fragment_match_store[cntr] = m
                    cntr += 1
                    if cntr % 100 == 0:
                        logger.info("Matched %d precursors", cntr)
                did_match_counter += counter
                combine_annotations(annotate_accumulator, annotater)
        else:
            logger.debug("Matching sequentially")
            for theoretical in theoretical_fragments:
                matches, counter, annotater = task_fn(theoretical)
                # results.extend(matches)
                for m in matches:
                    fragment_match_store[cntr] = m
                    cntr += 1
                did_match_counter += counter
                combine_annotations(annotate_accumulator, annotater)
    except Exception, e:
        logger.error(exc_info=e)
        raise e

    pool.close()
    pool.join()
    del pool
    if(cntr < 1):
        raise NoIonsMatchedException("No matches found from theoretical ions in MS2 deconvoluted results")
    logger.info("Saving.")
    fragment_match_store.commit()

    if(split_decon_data):
        try:
            logger.info("Splitting observed spectra between matched and unmatched")
            match, no_match = split_decon_data_by_index(data, did_match_counter)
            annotate_ions(match, annotate_accumulator)
            save_split_decon_data(match, no_match,
                                  "{root}.{tag}".format(root=splitext(decon_data)[0],
                                                        tag=tag))
        except Exception, e:
            logger.error(exc_info=e)

    return fragment_match_store.filename


def taskmain():
    import argparse
    logging.basicConfig(level='DEBUG')
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
