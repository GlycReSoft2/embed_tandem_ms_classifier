#!usr/bin/python

import json
import math
from os.path import splitext
from collections import Counter

import multiprocessing
import functools
import itertools

import logging
try:
    logging.basicConfig(stream=open("match-log.txt", 'w'), level=logging.INFO)
    pass
except:
    logging.basicConfig(level=logging.INFO)

from glycresoft_ms2_classification.error_code_interface import NoIonsMatchedException

from glycresoft_ms2_classification.data_io.bupid_topdown_deconvoluter import BUPIDYamlParser
from glycresoft_ms2_classification.data_io import default_loader, default_serializer, DeconIOBase

Proton = 1.007276035

ms1_tolerance_default = 10e-6
ms2_tolerance_default = 20e-6

MACHINE_ACQUISITION_RANGE = 3200 * 24


def save_interm(data, name):
    json.dump(data, open(name + ".json", 'wb'))


def Mergedicts(dicts_to_merge):
    # "firstrow" stores the data we will output. Items from other rows are incorporated
    # into it if they have a ppm error closer to 0, or if they have a row not existing in the first row.
    # in case the first dicts_to_merge is zero, we need a loop here. "allzero" is used
    # to check if all "rows in dicts_to_merge" are empty.
    # for dlist in dicts_to_merge:
    #     logging.debug("[")
    #     for d in dlist:
    #         logging.debug("\t{")
    #         for k, v in d.items():
    #             logging.debug("\t\t%s -> %s" % (k, v))
    #         logging.debug("\t}")
    #     logging.debug("]")
    allzero = True
    if len(dicts_to_merge) <= 1:
        if dicts_to_merge[0]:
            firstrow = []
            for items in dicts_to_merge[0]:
                sameMatchFound = False
                for items2 in firstrow:
                    if str(items['key']) == str(items2['key']):
                        sameMatchFound = True
                        if math.fabs(float(items['ppm_error'])) < math.fabs(float(items2['ppm_error'])):
                            items2 = items
                        break
                if sameMatchFound:
                    sameMatchFound = False
                else:
                    firstrow.append(items)
            return firstrow
        else:
            dicts_to_merge = []
            return dicts_to_merge
    for row in dicts_to_merge:
        if len(row) != 0:
            allzero = False
            firstrow = row
            break
    if allzero:
        dicts_to_merge = []
    # if the dicts_to_merge is empty anyway, just return the same
    # dicts_to_merge back.
        return dicts_to_merge
    else:
        for row in dicts_to_merge:
            if len(row) != 0:
                temprow = row
                firstrow = []
                for items in temprow:
                    sameMatchFound = False
                    for items2 in firstrow:
                        if str(items['key']) == str(items2['key']):
                            sameMatchFound = True
                            if math.fabs(float(items['ppm_error'])) < math.fabs(float(items2['ppm_error'])):
                                items2 = items
                            break
                    if sameMatchFound:
                        sameMatchFound = False
                    else:
                        firstrow.append(items)
                break
        # if the dicts_to_merge isn't empty, merge them.
        # Each of the "row in dicts_to_merge" should be something like this: [{'ppm_error':
        # -3.3899212497496274e-06, 'key': 'B3', 'obs_ion': 372.150116035},
        # {'ppm_error': -1.2204514502310785e-05, 'key': 'B4', 'obs_ion':
        # 532.175476035}, {'ppm_error': 1.812995934605501e-05, 'key': 'B4',
        # 'obs_ion': 532.191589035}, {'ppm_error': 1.4021538362178224e-07, 'key':
        # 'B5', 'obs_ion': 645.266113035}, {'ppm_error': 4.080114450361758e-06,
        # 'key': 'B7', 'obs_ion': 922.376038035}, {'ppm_error':
        # 1.0112467130843996e-05, 'key': 'B14', 'obs_ion': 1741.8025510349999}]
        for row in dicts_to_merge:
            # each of the "ReadingItems and firstrowReadingItems" should be
            # something like this: {'ppm_error': -3.3899212497496274e-06, 'key':
            # 'B3', 'obs_ion'}: 372.150116035}. If not, plz tell me.
            if len(row) != 0:
                for ReadingItems in row:
                    same_key_exists = False
                    for FR1 in range(len(firstrow)):
                        if str(ReadingItems['key']) == str(firstrow[FR1]['key']):
                            # if the same key exists, merge them.
                            same_key_exists = True
                            if math.fabs(float(ReadingItems['ppm_error'])) < math.fabs(float(firstrow[FR1]['ppm_error'])):
                                firstrow[FR1] = ReadingItems
                            break
                    # if the row doesn't exist in our final answer, add it.
                    if same_key_exists:
                        same_key_exists = False
                    else:
                        firstrow.append(ReadingItems)

    Finalfirstrow = []
    for items in firstrow:
        if items:
            Finalfirstrow.append(items)
    return Finalfirstrow

# TODO:
# Profile how scan_id spreads across merges
# Look at how much information can be discerned across time between two ambiguous matches
# Propagate this data along to downstream steps


def MergeRows(SourceData):
    # for row in range(len(SourceData)):
    #	print(SourceData[row]['MS1_Score'], SourceData[row]['bare_b_ions'], "\n")
    SourceData = sorted(SourceData, key=lambda k: k['Glycopeptide_identifier'])
    # MergedList is used to store our final result, the merged list.
    MergedList = []
    Storage = SourceData[0]
    newStorage = []
    for row in range(len(SourceData)):
        # newStorage is used to store the rows that have the same Glycopeptide identifier. Redeclaring it here to remove data from the last Glycopep ID.
        # input the first row.
        newStorage.append(Storage)
        #print(SourceData[row]['MS1_Score'], SourceData[row]['bare_b_ions'], "\n")
        if SourceData[row]['Glycopeptide_identifier'] == Storage['Glycopeptide_identifier']:
            logging.info(SourceData[row], Storage)
            Storage = SourceData[row]
            newStorage.append(Storage)
            # If this is the last row in SourceData, merge the rows now,
            # instead of doing it in the "else" statement.
            if (row + 1) == len(SourceData):
                # Declaring variables to store N P R T V X Y Z columns.
                oi = []
                bbi = []
                byi = []
                biwh = []
                yiwh = []
                bic = []
                yic = []
                si = []
                scan_id_range = set()
                # newrow is used to store our combined row
                newrow = []
                # input initial variables into newrow.
                newrow = newStorage[0]

                for row2 in newStorage:
                    oi.append(row2["Oxonium_ions"])
                    bbi.append(row2["bare_b_ions"])
                    byi.append(row2["bare_y_ions"])
                    biwh.append(row2["b_ions_with_HexNAc"])
                    yiwh.append(row2["y_ions_with_HexNAc"])
                    bic.append(row2["b_ion_coverage"])
                    yic.append(row2["y_ion_coverage"])
                    si.append(row2["Stub_ions"])
                    scan_id_range.add(row2["scan_id"])
                # Change the N P R T V X Y Z columns' values in newrow.
                #print(newrow['MS1_Score'], yic, "\n")
                newrow["Oxonium_ions"] = Mergedicts(oi)
                newrow["bare_b_ions"] = Mergedicts(bbi)
                newrow["bare_y_ions"] = Mergedicts(byi)
                newrow["b_ions_with_HexNAc"] = Mergedicts(biwh)
                newrow["y_ions_with_HexNAc"] = Mergedicts(yiwh)
                newrow["b_ion_coverage"] = Mergedicts(bic)
                newrow["y_ion_coverage"] = Mergedicts(yic)
                newrow["Stub_ions"] = Mergedicts(si)
                newrow["scan_id_range"] = list(scan_id_range)
                # Now put the newrow into the final output
                #print (newrow['MS1_Score'], newrow["y_ion_coverage"], "\n\n")
                MergedList.append(newrow)
                newStorage = []
                break

        else:
            # The new row has another glycopep identifier, so we are storing it
            # in Storage for the next loop.
            Storage = SourceData[row]
            # Lets start merging the rows now.
            # Declaring variables to store N P R T V X Y Z columns.
            oi = []
            bbi = []
            byi = []
            biwh = []
            yiwh = []
            bic = []
            yic = []
            si = []
            scan_id_range = set()
            # newrow is used to store our combined row
            newrow = []
            # input initial variables into newrow.
            newrow = newStorage[0]
            for row2 in newStorage:
                oi.append(row2["Oxonium_ions"])
                bbi.append(row2["bare_b_ions"])
                byi.append(row2["bare_y_ions"])
                biwh.append(row2["b_ions_with_HexNAc"])
                yiwh.append(row2["y_ions_with_HexNAc"])
                bic.append(row2["b_ion_coverage"])
                yic.append(row2["y_ion_coverage"])
                si.append(row2["Stub_ions"])
                scan_id_range.add(row2["scan_id"])
            # Change the N P R T V X Y Z columns' values in newrow.
            #print(newrow['MS1_Score'], yic, "\n")
            newrow["Oxonium_ions"] = Mergedicts(oi)
            newrow["bare_b_ions"] = Mergedicts(bbi)
            newrow["bare_y_ions"] = Mergedicts(byi)
            newrow["b_ions_with_HexNAc"] = Mergedicts(biwh)
            newrow["y_ions_with_HexNAc"] = Mergedicts(yiwh)
            newrow["b_ion_coverage"] = Mergedicts(bic)
            newrow["y_ion_coverage"] = Mergedicts(yic)
            newrow["Stub_ions"] = Mergedicts(si)
            newrow["scan_id_range"] = list(scan_id_range)
            #print (newrow['MS1_Score'], newrow["y_ion_coverage"], "\n\n")
            # Now put the newrow into the final output
            newStorage = []
            MergedList.append(newrow)
            continue
    return MergedList


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
                            {"ion": (ob_ions.mass + Proton), "ppm_error": oxonium_ppm * 1e6, "key": ox_ion})

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
                             "ppm_error": tandem_ppm * 1e6})
                        all_b_ions.append(
                            {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                             "ppm_error": tandem_ppm * 1e6})

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
                             "ppm_error": tandem_ppm * 1e6})
                        all_b_ions.append(
                            {"obs_ion": (obs_ions.mass + Proton),
                             "key": theo_ions["key"].split("+")[0],
                             "ppm_error": tandem_ppm * 1e6})

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
                             "ppm_error": tandem_ppm * 1e6})
                        all_y_ions.append(
                            {"obs_ion": (obs_ions.mass + Proton), "key": theo_ions["key"],
                             "ppm_error": tandem_ppm * 1e6})

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
                             "ppm_error": tandem_ppm * 1e6})
                        all_y_ions.append(
                            {"obs_ion": (obs_ions.mass + Proton),
                             "key": theo_ions["key"].split("+")[0],
                             "ppm_error": tandem_ppm * 1e6})

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
                             "ppm_error": tandem_ppm * 1e6})
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
    grouped = itertools.groupby(matches, lambda x: x["Glycopeptide_identifier"])
    merged_matches = []
    for glycopeptide, matches in grouped:
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
    best_matches = []
    fabs = math.fabs
    grouped_by_key = itertools.groupby(itertools.chain.from_iterable(matches), lambda x: x["key"])
    for key, matched_key in grouped_by_key:
        best_match = matched_key.next()
        best_ppm = fabs(best_match["ppm_error"])
        for match in matched_key:
            if fabs(match["ppm_error"]) < best_ppm:
                best_match = match
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
    if outfile is None:
        outfile = splitext(splitext(theo_fragment_file)[0])[0] + '.match_frags'

    pool = multiprocessing.Pool(n_processes)

    data = decon_format_lookup[decon_data_format](decon_data)
    print("Loading observed ions complete")

    if isinstance(theo_fragment_file, (str, file)):
        theoretical_search_space = json.load(open(theo_fragment_file))
        if outfile is None:
            outfile = splitext(splitext(theo_fragment_file)[0])[0] + '.match_frags'
    else:
        theoretical_search_space = theo_fragment_file

    metadata = theoretical_search_space["metadata"]
    metadata["ms1_ppm_tolerance"] = ms1_tolerance
    metadata["ms2_ppm_tolerance"] = ms2_tolerance

    if outfile is None:
        outfile = splitext(metadata["ms1_output_file"])[0] + '.match_frags'

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
    task_fn = functools.partial(match_observed_to_theoretical, ms1_tolerance=ms1_tolerance,
                                ms2_tolerance=ms2_tolerance, observed_ions=data)
    matcheing_process = pool.map(task_fn, theoretical_fragments, chunksize=len(theoretical_fragments)/n_processes)
    for matches, counter in matcheing_process:
        results.extend(matches)
        did_match_counter += counter

    if(len(results) < 1):
        raise NoIonsMatchedException(
            "No matches found from theoretical ions in MS2 deconvoluted results")

    save_unmerged = multiprocessing.Process(target=save_interm, args=(results, "unmerged-matches"))
    save_unmerged.start()

    merged_results = merge_by_glycopeptide_sequence(results)

    results = {
        "metadata": metadata,
        "matched_ions": merged_results
    }
    print("Matching done")
    f = open(outfile + '.json', 'wb')
    json.dump(results, f)
    f.close()

    if(split_decon_data):
        print("Splitting")
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
