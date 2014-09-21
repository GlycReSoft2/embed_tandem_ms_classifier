#!usr/bin/python

import ast
import csv
#import json
import math
#import os
from os.path import splitext
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

ms1_tolerance_default = 0.00001
ms2_tolerance_default = 0.00002


def Mergedicts(dicts_to_merge):
    # "firstrow" stores the data we will output. Items from other rows are incorporated
    # into it if they have a ppm error closer to 0, or if they have a row not existing in the first row.
    # in case the first dicts_to_merge is zero, we need a loop here. "allzero" is used
    # to check if all "rows in dicts_to_merge" are empty.
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
    # if the dicts_to_merge is empty anyway, just return the same dicts_to_merge back.
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
    isB = 0
    for items in firstrow:
        if items:
            Finalfirstrow.append(items)
    return Finalfirstrow


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
            #print (newrow['MS1_Score'], newrow["y_ion_coverage"], "\n\n")
            # Now put the newrow into the final output
            newStorage = []
            MergedList.append(newrow)
            continue
    return MergedList


def match_frags(theo_fragment_file, yaml_data, ms1_tolerance=ms1_tolerance_default,
                ms2_tolerance=ms2_tolerance_default, outfile=None):
    '''
    :param theo_fragment_file: path to file containing all theoretical peptides
    :param yaml_data: path to yaml file from earlier MS1 analysis
    :param outfile: path to the file to write output to. Defaults to theo_fragment_file + ".match_frags".
    '''
    if outfile is None:
        outfile = splitext(splitext(theo_fragment_file)[0])[0] + '.match_frags'

    #filename = "Transferrin-ions2.csv"
    f_csv = csv.DictReader(open(theo_fragment_file))
    #datafile = "KK_20130630_007_group_ms2.yaml"
    stream = open(yaml_data, 'r')
    data = []
    data.append(yaml.load(stream))

    results = []

    # reading each line in theoretical ions file generated from gly1 results
    for lines in f_csv:

        for rows in data:  # reading yaml data scan by scan

            for more in rows['peaks']:

                real_ions = []
                # all ions in scan assigned to list real_ions
                real_ions = more['mass']
                scan_id = more['scans'][0]['id']
                #print ("intensity array", more['intensity'])
                #print ("number of elements" , more['num'])

                # this takes the information from the first scan in case there
                # are multiple merged scans
                info = more['scans'][:1]
                for num in info:

                    #print ("mz", info['mz'])
                    charge = num['z']
                    mass = num['mz']
                    ntr_mass = ((mass * charge) - (charge * Proton))
                    #print ("ntr-mass", ntr_mass)

                    obs_mass = float(lines['Obs_Mass'])
                    ppm = ((ntr_mass - obs_mass) / ntr_mass)

                if math.fabs(ppm) <= ms1_tolerance:
                    Oxonium_ions = ast.literal_eval(lines['Oxonium_ions'])

                    oxoniums = []
                    for ox_ion in Oxonium_ions:
                        for ob_ions in real_ions:
                            oxonium_ppm = (
                                ((ob_ions + Proton) - ox_ion) / (ob_ions + Proton))
                            if math.fabs(oxonium_ppm) <= ms2_tolerance:
                                oxoniums.append(
                                    {"ion": (ob_ions + Proton), "ppm_error": oxonium_ppm, "key": ox_ion})
                                #print ("\n", ob_ions, "ppm error:", oxonium_ppm)
                            else:
                                pass
                    #oxo_found = len(oxoniums)

                    # checking for b and y ions and ions with HexNAc:
                    b_ions = ast.literal_eval(lines['bare_b_ions'])
                    b_len = float(len(b_ions))
                    # print "\n, b ions:"
                    b_type = []
                    all_b_ions = []
                    for theo_ions in b_ions:
                        frag_mass = float(theo_ions["mass"])
                        # change this to proton once the error is found in
                        # Han's code
                        prot_ion = frag_mass - Proton
                        for obs_ions in real_ions:
                            tandem_ppm = float(
                                (obs_ions - (prot_ion)) / obs_ions)
                            if math.fabs(tandem_ppm) <= ms2_tolerance:
                                b_type.append(
                                    {"obs_ion": (obs_ions + Proton), "key": theo_ions["key"], "ppm_error": tandem_ppm})
                                all_b_ions.append(
                                    {"obs_ion": (obs_ions + Proton), "key": theo_ions["key"], "ppm_error": tandem_ppm})
                                #print (obs_ions, theo_ions["mass"], theo_ions["key"], tandem_ppm)
                            else:
                                pass

                    #len_found = float(len(b_type))
                    # if len_found != 0:
                    #	b_found = (len_found/b_len)*100
                    # else:
                    #	b_found = 0

                    b_HexNAc_ions = ast.literal_eval(
                        lines['b_ions_with_HexNAc'])
                    b_HexNAc_len = float(len(b_HexNAc_ions))
                    # print "\n, b ions:"
                    b_HexNAc_type = []
                    for theo_ions in b_HexNAc_ions:
                        frag_mass = float(theo_ions["mass"])
                        # change this to proton once the error is found in
                        # Han's code
                        prot_ion = frag_mass - Proton
                        for obs_ions in real_ions:
                            tandem_ppm = float(
                                (obs_ions - (prot_ion)) / obs_ions)
                            if math.fabs(tandem_ppm) <= ms2_tolerance:
                                b_HexNAc_type.append(
                                    {"obs_ion": (obs_ions + Proton), "key": theo_ions["key"], "ppm_error": tandem_ppm})
                                all_b_ions.append(
                                    {"obs_ion": (obs_ions + Proton), "key": theo_ions["key"].split("+")[0], "ppm_error": tandem_ppm})
                                #print (obs_ions, theo_ions["mass"], theo_ions["key"], tandem_ppm)
                            else:
                                pass
                    #b_HexNAc_len_found = len(b_HexNAc_type)

                    y_ions = list(ast.literal_eval(lines['bare_y_ions']))
                    # print "\n, y ions:"
                    y_len = float(len(y_ions))
                    y_type = []
                    all_y_ions = []
                    for theo_ions in y_ions:
                        frag_mass = float(theo_ions["mass"])
                        # change this to proton once the error is found in
                        # Han's code
                        prot_ion = frag_mass - Proton
                        for obs_ions in real_ions:
                            tandem_ppm = float(
                                (obs_ions - (prot_ion)) / obs_ions)
                            if math.fabs(tandem_ppm) <= ms2_tolerance:
                                y_type.append(
                                    {"obs_ion": (obs_ions + Proton), "key": theo_ions["key"], "ppm_error": tandem_ppm})
                                all_y_ions.append(
                                    {"obs_ion": (obs_ions + Proton), "key": theo_ions["key"], "ppm_error": tandem_ppm})
                                #print (obs_ions, theo_ions["mass"], theo_ions["key"], tandem_ppm, "\n")
                            else:
                                pass
                    #len_found = float(len(y_type))
                    # if len_found != 0:
                        #y_found = (len_found/y_len)*100
                    # elif len_found ==0:
                        #y_found = 0

                    y_HexNAc_ions = ast.literal_eval(
                        lines['y_ions_with_HexNAc'])
                    y_HexNAc_len = float(len(y_HexNAc_ions))
                    # print "\n, b ions:"
                    y_HexNAc_type = []
                    for theo_ions in y_HexNAc_ions:
                        frag_mass = float(theo_ions["mass"])
                        # change this to proton once the error is found in
                        # Han's code
                        prot_ion = frag_mass - Proton
                        for obs_ions in real_ions:
                            tandem_ppm = float(
                                (obs_ions - (prot_ion)) / obs_ions)
                            if math.fabs(tandem_ppm) <= ms2_tolerance:
                                y_HexNAc_type.append(
                                    {"obs_ion": (obs_ions + Proton), "key": theo_ions["key"], "ppm_error": tandem_ppm})
                                all_y_ions.append(
                                    {"obs_ion": (obs_ions + Proton), "key": theo_ions["key"].split("+")[0], "ppm_error": tandem_ppm})
                                #print (obs_ions, theo_ions["mass"], theo_ions["key"], tandem_ppm)
                            else:
                                pass
                    #y_HexNAc_len_found = len(y_HexNAc_type)

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
                                    {"obs_ion": (obs_ions + Proton), "key": theo_ions["key"], "ppm_error": tandem_ppm})
                                #print (obs_ions, theo_ions["mass"], theo_ions["key"], tandem_ppm, "\n")
                            else:
                                pass
                    # if len(stub_type) != 0:
                    #	stub_found = len (stub_type)
                    # else:
                    #	stub_found = 0
                    results.append(
                        {"MS1_Score": lines["MS1_Score"], "Obs_Mass": lines["Obs_Mass"], "Calc_mass": lines["Calc_mass"],
                         "ppm_error": lines["ppm_error"], "Peptide": lines["Peptide"], "Peptide_mod": lines["Peptide_mod"],
                         "Glycan": lines["Glycan"],
                         "vol": lines["vol"], "glyco_sites": lines["glyco_sites"], "startAA": lines["startAA"], "endAA": lines["endAA"],
                         "Seq_with_mod": lines["Seq_with_mod"], "Glycopeptide_identifier": lines["Seq_with_mod"] + lines["Glycan"],
                         "Oxonium_ions": oxoniums, "bare_b_ions": b_type, "possible_b_ions_HexNAc": len(b_HexNAc_ions),
                         "total_b_ions": b_len, "bare_y_ions": y_type, "possible_y_ions_HexNAc": len(y_HexNAc_ions),
                         "total_y_ions": y_len, "b_ions_with_HexNAc": b_HexNAc_type, "y_ions_with_HexNAc": y_HexNAc_type,
                         "b_ion_coverage": all_b_ions, "y_ion_coverage": all_y_ions, "Stub_ions": stub_type,
                         "scan_id": scan_id
                         })
                else:
                    pass

    if(len(results) < 1):
        raise NoIonsMatchedException("No matches found from theoretical ions in MS2 deconvoluted results")

    merged_results = MergeRows(results)
    # for i in merged_results:
    #     if i["Glycopeptide_identifier"] == "DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;5;4;2;0]":
    # print ("l3", i['b_ion_coverage'])
    #         pass

    keys = [
        "MS1_Score", "Obs_Mass", "Calc_mass", "ppm_error", "Peptide", "Peptide_mod", "Glycan", "vol", "glyco_sites",
        "startAA", "endAA", "Seq_with_mod", "Glycopeptide_identifier", "Oxonium_ions",
        "bare_b_ions", "possible_b_ions_HexNAc", "total_b_ions", "bare_y_ions", "possible_y_ions_HexNAc",
        "total_y_ions", "b_ions_with_HexNAc", "y_ions_with_HexNAc", "b_ion_coverage", "y_ion_coverage",
        "Stub_ions", "scan_id"]
    f = open(outfile + '.csv', 'wb')
    dict_writer = csv.DictWriter(f, keys)
    dict_writer.writer.writerow(keys)
    dict_writer.writerows(merged_results)
    f.close()
    return outfile + '.csv'




if __name__ == '__main__':
    import sys
    match_frags(sys.argv[1], sys.argv[2])
