import math
import multiprocessing
import functools
import itertools
import logging

from os.path import splitext, basename
from collections import Counter, defaultdict

import sqlitedict

from glycresoft_ms2_classification.structure.composition import Composition
from glycresoft_ms2_classification.error_code_interface import NoIonsMatchedException
from glycresoft_ms2_classification.utils import collectiontools
from glycresoft_ms2_classification.ms.bupid_topdown_deconvoluter import BUPIDYamlParser
from glycresoft_ms2_classification.ms import default_loader, default_serializer
from glycresoft_ms2_classification.ms import DeconIOBase, MSMSSqlDB, ObservedPrecursorSpectrum
from glycresoft_ms2_classification.utils.parallel_opener import ParallelParser

from .results_merger import Merger

logger = logging.getLogger(__name__)

use_cython = True
try:
    from glycresoft_ms2_classification.ms import ion_matching
    c_match_observed_to_theoretical_sql = ion_matching.match_observed_wrap
except:
    use_cython = False

Proton = Composition("H+").mass
ms1_tolerance_default = 10e-6
ms2_tolerance_default = 20e-6


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
                         "key": kind, "scan_id": tandem_ms_ind})
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
                         "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
                    all_b_ions.append(
                        {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                         "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
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
                         "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
                    all_b_ions.append(
                        {"obs_ion": (obs_ions.mass + Proton),
                         "key": theo_ions["key"].split("+")[0],
                         "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
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
                         "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
                    all_y_ions.append(
                        {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                         "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
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
                         "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
                    all_y_ions.append(
                        {"obs_ion": (obs_ions.mass + Proton),
                         "key": theo_ions["key"].split("+")[0],
                         "ppm_error": tandem_ppm * 1e6, "scan_id": tandem_ms_ind})
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


def dispatch_task(theoretical, observed_ions_conn_string, ms1_tolerance, ms2_tolerance,
                  worker_fn, out_queue):
    res = worker_fn(theoretical, observed_ions_conn_string, ms1_tolerance, ms2_tolerance)
    out_queue.put(res)
    return len(res[0])


class MatchDispatcher(object):
    def __init__(self, fragment_db_file, obs_data_db_file, n_processes=4,
                 ms1_tolerance=ms1_tolerance_default,
                 ms2_tolerance=ms2_tolerance_default, chunk_size=500):
        self.fragment_db_file = fragment_db_file
        self.obs_data_db_file = obs_data_db_file
        self.n_processes = n_processes
        self.ms1_tolerance = ms1_tolerance
        self.ms2_tolerance = ms2_tolerance
        self.chunk_size = chunk_size
        self.worker_pool = multiprocessing.Pool(self.n_processes)
        merge_file_base = fragment_db_file + "_" + basename(splitext(obs_data_db_file)[0])
        self.results_merger = Merger(merge_file_base)
        match_fn = match_observed_to_theoretical_sql
        if use_cython:
            logger.info("Using Cython")
            match_fn = c_match_observed_to_theoretical_sql
        self.task_fn = functools.partial(dispatch_task, worker_fn=match_fn,
                                         out_queue=self.results_merger.data_queue,
                                         ms1_tolerance=self.ms1_tolerance,
                                         ms2_tolerance=self.ms2_tolerance,
                                         observed_ions_conn_string=obs_data_db_file)

    def run(self):
        self.results_merger.start_worker()
        logger.debug("Matching concurrently")
        theoretical_search_space = sqlitedict.SqliteDict(
            self.fragment_db_file, tablename="theoretical_search_space")
        theoretical_fragments = theoretical_search_space.itervalues()

        matching_process = self.worker_pool.imap_unordered(
            self.task_fn, theoretical_fragments, chunksize=self.chunk_size)
        for match_count in matching_process:
            pass

        self.results_merger.stop_worker()
        logger.info("Matching Complete")
