#!usr/bin/python

import ast
import csv
import math
from os.path import splitext
from collections import Counter
from itertools import izip


import logging
try:
    logging.basicConfig(stream=open("match-log.txt", 'w'), level=logging.INFO)
except:
    logging.basicConfig(level=logging.INFO)
import yaml

from error_code_interface import NoIonsMatchedException

Carbon = 12.00000
Hydrogen = 1.0078250350
Nitrogen = 14.0030740000
Oxygen = 15.9949146300
Sulfur = 31.9720707000
Phosphorus = 30.9737620000
Iodine = 126.904473
Chlorine = 34.9688527300
Fluorine = 18.9984032200
Bromine = 78.918361000
Sodium = 22.989767700
Potassium = 38.9637069000
Calcium = 39.9625912000
Electron = 0.000549
Water = 18.0105647000
Proton = 1.007276035

ms1_tolerance_default = 10e-6
ms2_tolerance_default = 20e-6

MACHINE_ACQUISITION_RANGE = 3200 * 24


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
    matched = []
    no_matched = []
    for i, scan in enumerate(decon_data['peaks']):
        count = splitting_index[i]
        scan["__matches"] = count
        if count != 0:
            matched.append(scan)
        else:
            no_matched.append(scan)
    return matched, no_matched


def match_frags(theo_fragment_file, decon_data, ms1_tolerance=ms1_tolerance_default,
                ms2_tolerance=ms2_tolerance_default, split_decon_data=False, outfile=None):
    '''
    :param theo_fragment_file: path to file containing all theoretical peptides
    :param decon_data: path to yaml file from earlier MS1 analysis
    :param ms1_tolerance: float tolerance for matching intact mass
    :param ms2_tolerance: float tolerance for matching fragment ions
    :param split_decon_data: bool flag whether to emit a file for all scans without matches
                             and a file for all scans with matches.
    :param outfile: path to the file to write output to. Defaults to theo_fragment_file + ".match_frags".
    '''
    if outfile is None:
        outfile = splitext(splitext(theo_fragment_file)[0])[0] + '.match_frags'

    f_csv = csv.DictReader(open(theo_fragment_file))
    stream = open(decon_data, 'r')
    data = []
    try:
        loader = yaml.CLoader(stream)
    except:
        loader = yaml.Loader(stream)
    data = (loader.get_data())

    results = []
    did_match_cache = Counter()
    # reading each line in theoretical ions file generated from gly1 results
    for line_no, lines in enumerate(f_csv):
        try:
            obs_mass = float(lines['Calc_mass'])
        except:
            print(lines["Calc_mass"], line_no)
            print(lines)
            raise
        for tandem_ms_ind, more in enumerate(data['peaks']):
            more['_iterkey'] = tandem_ms_ind
            real_ions = []
            # all ions in scan assigned to list real_ions
            # real_ions = more['mass']
            scan_id = more['scans'][0]['id']

            # this takes the information from the first scan in case there
            # are multiple merged scans
            num = more['scans'][0]
            # for num in info:
            charge = num['z']
            mass = num['mz']
            ntr_mass = ((mass * charge) - (charge * Proton))
            precursor_ppm = ((ntr_mass - obs_mass) / ntr_mass)
            
            if math.fabs(precursor_ppm) <= ms1_tolerance:
                did_match_cache[tandem_ms_ind] += 1
                Oxonium_ions = ast.literal_eval(lines['Oxonium_ions'])

                # Try to limit available ions to match by the precursor charge
                try:
                    real_ions = [ion[0] for ion
                                 in izip(more["mass"], more["z"], more["intensity"])
                                 if ion[1] <= charge and ion[0] < MACHINE_ACQUISITION_RANGE]
                except Exception, e:
                    logging.info(
                        "An error occurred while trying to filter ions")
                    logging.exception(e)
                    real_ions = more["mass"]

                oxoniums = []
                for ox_ion in Oxonium_ions:
                    for ob_ions in real_ions:
                        oxonium_ppm = (
                            ((ob_ions + Proton) - ox_ion) / (ob_ions + Proton))
                        if math.fabs(oxonium_ppm) <= ms2_tolerance:
                            oxoniums.append(
                                {"ion": (ob_ions + Proton), "ppm_error": oxonium_ppm * 1e6, "key": ox_ion})

                # checking for b and y ions and ions with HexNAc:
                b_ions = ast.literal_eval(lines['bare_b_ions'])
                b_len = float(len(b_ions))
                b_type = []
                all_b_ions = []
                for theo_ions in b_ions:
                    frag_mass = float(theo_ions["mass"])
                    prot_ion = frag_mass - Proton
                    for obs_ions in real_ions:
                        tandem_ppm = float(
                            (obs_ions - (prot_ion)) / obs_ions)
                        if math.fabs(tandem_ppm) <= ms2_tolerance:
                            b_type.append(
                                {"obs_ion": (obs_ions + Proton), "key": theo_ions["key"],
                                 "ppm_error": tandem_ppm * 1e6})
                            all_b_ions.append(
                                {"obs_ion": (obs_ions + Proton), "key": theo_ions["key"],
                                 "ppm_error": tandem_ppm * 1e6})

                b_HexNAc_ions = ast.literal_eval(
                    lines['b_ions_with_HexNAc'])
                b_HexNAc_type = []
                for theo_ions in b_HexNAc_ions:
                    frag_mass = float(theo_ions["mass"])
                    prot_ion = frag_mass - Proton
                    for obs_ions in real_ions:
                        tandem_ppm = float(
                            (obs_ions - (prot_ion)) / obs_ions)
                        if math.fabs(tandem_ppm) <= ms2_tolerance:
                            b_HexNAc_type.append(
                                {"obs_ion": (obs_ions + Proton), "key": theo_ions["key"],
                                 "ppm_error": tandem_ppm * 1e6})
                            all_b_ions.append(
                                {"obs_ion": (obs_ions + Proton),
                                 "key": theo_ions["key"].split("+")[0],
                                 "ppm_error": tandem_ppm * 1e6})

                y_ions = list(ast.literal_eval(lines['bare_y_ions']))
                y_len = float(len(y_ions))
                y_type = []
                all_y_ions = []
                for theo_ions in y_ions:
                    frag_mass = float(theo_ions["mass"])
                    prot_ion = frag_mass - Proton
                    for obs_ions in real_ions:
                        tandem_ppm = float(
                            (obs_ions - (prot_ion)) / obs_ions)
                        if math.fabs(tandem_ppm) <= ms2_tolerance:
                            y_type.append(
                                {"obs_ion": (obs_ions + Proton), "key": theo_ions["key"],
                                 "ppm_error": tandem_ppm * 1e6})
                            all_y_ions.append(
                                {"obs_ion": (obs_ions + Proton), "key": theo_ions["key"],
                                 "ppm_error": tandem_ppm * 1e6})

                y_HexNAc_ions = ast.literal_eval(
                    lines['y_ions_with_HexNAc'])
                y_HexNAc_type = []
                for theo_ions in y_HexNAc_ions:
                    frag_mass = float(theo_ions["mass"])
                    prot_ion = frag_mass - Proton
                    for obs_ions in real_ions:
                        tandem_ppm = float(
                            (obs_ions - (prot_ion)) / obs_ions)
                        if math.fabs(tandem_ppm) <= ms2_tolerance:
                            y_HexNAc_type.append(
                                {"obs_ion": (obs_ions + Proton), "key": theo_ions["key"],
                                 "ppm_error": tandem_ppm * 1e6})
                            all_y_ions.append(
                                {"obs_ion": (obs_ions + Proton),
                                 "key": theo_ions["key"].split("+")[0],
                                 "ppm_error": tandem_ppm * 1e6})

                # checking for stub ions
                stub_ions = ast.literal_eval(lines['pep_stub_ions'])
                stub_type = []
                for theo_ions in stub_ions:
                    frag_mass = float(theo_ions["mass"])
                    prot_ion = frag_mass - Proton
                    for obs_ions in real_ions:
                        tandem_ppm = float(
                            (obs_ions - (prot_ion)) / obs_ions)
                        if math.fabs(tandem_ppm) <= ms2_tolerance:
                            stub_type.append(
                                {"obs_ion": (obs_ions + Proton), "key": theo_ions["key"],
                                 "ppm_error": tandem_ppm * 1e6})

                results.append(
                    {"MS1_Score": lines["MS1_Score"], "Obs_Mass": lines["Obs_Mass"],
                     "Calc_mass": lines["Calc_mass"],
                     "ppm_error": lines["ppm_error"], "Peptide": lines["Peptide"],
                     "Peptide_mod": lines["Peptide_mod"],
                     "Glycan": lines["Glycan"],
                     "vol": lines["vol"], "glyco_sites": lines["glyco_sites"],
                     "startAA": lines["startAA"], "endAA": lines["endAA"],
                     "Seq_with_mod": lines["Seq_with_mod"],
                     "Glycopeptide_identifier": lines["Seq_with_mod"] + lines["Glycan"],
                     "Oxonium_ions": oxoniums, "bare_b_ions": b_type, "possible_b_ions_HexNAc": len(b_HexNAc_ions),
                     "total_b_ions": b_len, "bare_y_ions": y_type, "possible_y_ions_HexNAc": len(y_HexNAc_ions),
                     "total_y_ions": y_len, "b_ions_with_HexNAc": b_HexNAc_type,
                     "y_ions_with_HexNAc": y_HexNAc_type,
                     "b_ion_coverage": all_b_ions, "y_ion_coverage": all_y_ions, "Stub_ions": stub_type,
                     "scan_id": scan_id,
                     "scan_id_range": []
                     })
            else:
                pass

    if(len(results) < 1):
        raise NoIonsMatchedException(
            "No matches found from theoretical ions in MS2 deconvoluted results")

    merged_results = MergeRows(results)

    keys = [
        "MS1_Score", "Obs_Mass", "Calc_mass", "ppm_error", "_old_ppm_error", "Peptide",
        "Peptide_mod", "Glycan", "vol", "glyco_sites",
        "startAA", "endAA", "Seq_with_mod", "Glycopeptide_identifier", "Oxonium_ions",
        "bare_b_ions", "possible_b_ions_HexNAc", "total_b_ions", "bare_y_ions", "possible_y_ions_HexNAc",
        "total_y_ions", "b_ions_with_HexNAc", "y_ions_with_HexNAc", "b_ion_coverage", "y_ion_coverage",
        "Stub_ions", "scan_id", "scan_id_range"
    ]
    f = open(outfile + '.csv', 'wb')
    dict_writer = csv.DictWriter(f, keys)
    dict_writer.writer.writerow(keys)
    dict_writer.writerows(merged_results)
    f.close()

    if(split_decon_data):
        print("Splitting decon_data by matching criterion")
        match, no_match = split_decon_data_by_index(data, did_match_cache)
        decon_path = splitext(decon_data)[0]
        match_fh = open(decon_path + '-matches.yaml', 'wb')
        try:
            dumper = yaml.CDumper(match_fh)
        except:
            dumper = yaml.Dumper(match_fh)
        dumper.open()
        dumper.represent(match)
        match_fh.close()
        no_match_fh = open(decon_path + '-no-matches.yaml', 'wb')
        try:
            dumper = yaml.CDumper(no_match_fh)
        except:
            dumper = yaml.Dumper(no_match_fh)
        dumper.open()
        dumper.represent(no_match)

    return outfile + '.csv'


if __name__ == '__main__':
    import sys
    match_frags(sys.argv[1], sys.argv[2])
