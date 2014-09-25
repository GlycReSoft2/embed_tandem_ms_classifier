import ast
import csv
import json
#import math
from os.path import splitext
import re


def main(matched_ions_file, output_file=None):
    '''
    Scores observed fragment to theoretical ion pairings, eliminates ambiguity using tandem MS2 data
    :param matched_ions_file: Path to file created by match_ions.py, matching theoretical ions to observed spectra
    :param output_file: Path to file to write results to. Defaults to `matched_ions_file` + ".processed"
    '''
    if output_file is None:
        output_file = splitext(
            splitext(matched_ions_file)[0])[0] + ".processed"
    results = csv.DictReader(open(matched_ions_file))
    scored_results = []

    for rows in results:

        MS2_score = 0.0

        total_bare_b = float(rows["total_b_ions"])
        # print total_bare_b

        total_bare_y = float(rows["total_y_ions"])

        all_b_found = ast.literal_eval(rows["b_ion_coverage"])
        num_b_found = float(len(all_b_found))
        percent_b = (num_b_found / total_bare_b)

        all_y_found = ast.literal_eval(rows["y_ion_coverage"])
        num_y_found = float(len(all_y_found))
        percent_y = (num_y_found / total_bare_y)

        b_HexNAc = ast.literal_eval(rows["b_ions_with_HexNAc"])
        b_HexNAc_found = len(b_HexNAc)
        total_possible_b_HexNAc = float(
            ast.literal_eval(rows["possible_b_ions_HexNAc"]))
        if total_possible_b_HexNAc == 0:
            percent_b_HexNAc = 0
        else:
            percent_b_HexNAc = b_HexNAc_found / total_possible_b_HexNAc

        y_HexNAc = ast.literal_eval(rows["y_ions_with_HexNAc"])
        y_HexNAc_found = len(y_HexNAc)
        total_possible_y_HexNAc = float(
            ast.literal_eval(rows["possible_y_ions_HexNAc"]))
        if total_possible_y_HexNAc == 0:
            percent_y_HexNAc = 0
        else:
            percent_y_HexNAc = y_HexNAc_found / total_possible_y_HexNAc

        stubs = ast.literal_eval(rows["Stub_ions"])
        stubs_found = len(stubs)

        Oxoniums = ast.literal_eval(rows["Oxonium_ions"])
        oxoniums_found = len(Oxoniums)

        glycans = re.findall('\\b\\d+\\b', rows["Glycan"])
        NeuAc = int(glycans[3])

        # checking for oxonium ions, if less than 2, MS2 score = 0
        if len(Oxoniums) < 2:
            MS2_score = 0.0
        else:                   # proceed only if 2 or more oxoniums are found
            MS2_score += (0.02 * oxoniums_found)

            MS2_score += (
                (percent_b * 2.00) + (percent_y * 2.00) + (0.01 * stubs_found))

            # decrease the score if NeuAc oxoniums are seen in a Non-NeuAc
            # composition
            if NeuAc != 0:
                pass

            else:
                if any(d['key'] == 274.092 or d["key"] == 292.1026 for d in Oxoniums):
                    MS2_score -= 0.2

                else:
                    pass

        if (y_HexNAc_found == 0 and b_HexNAc_found == 0):
            pass
        else:
            MS2_score += (
                (percent_b_HexNAc * 0.02) + (percent_y_HexNAc * 0.02))

        scored_results.append({"MS1_Score": rows["MS1_Score"],
                               "Obs_Mass": rows["Obs_Mass"],
                               "Calc_mass": rows["Calc_mass"],
                               "ppm_error": rows["ppm_error"],
                               "Peptide": rows["Peptide"],
                               "Peptide_mod": rows["Peptide_mod"],
                               "Glycan": rows["Glycan"],
                               "vol": rows["vol"],
                               "glyco_sites": rows["glyco_sites"],
                               "startAA": rows["startAA"],
                               "endAA": rows["endAA"],
                               "Seq_with_mod": rows["Seq_with_mod"],
                               "Glycopeptide_identifier": rows["Seq_with_mod"] + rows["Glycan"],
                               "Oxonium_ions": json.dumps(Oxoniums),
                               "bare_b_ions": json.dumps(ast.literal_eval(rows["bare_b_ions"])),
                               "total_b_ions_possible": total_bare_b,
                               "bare_y_ions": json.dumps(ast.literal_eval(rows["bare_y_ions"])),
                               "total_y_ions_possible": total_bare_y,
                               "b_ions_with_HexNAc": json.dumps(b_HexNAc),
                               "y_ions_with_HexNAc": json.dumps(y_HexNAc),
                               "b_ion_coverage": json.dumps(all_b_found),
                               "y_ion_coverage": json.dumps(all_y_found),
                               "Stub_ions": json.dumps(stubs),
                               "percent_b_ion_with_HexNAc_coverage": percent_b_HexNAc,
                               "percent_y_ion_with_HexNAc_coverage": percent_y_HexNAc,
                               "percent_b_ion_coverage": percent_b,
                               "percent_y_ion_coverage": percent_y,
                               "#_of_stubs_found": stubs_found,
                               "scan_id": rows['scan_id']
                               })

    keys = [
        "MS1_Score",
        "Obs_Mass",
        "Calc_mass",
        "ppm_error",
        "Peptide",
        "Peptide_mod",
        "Glycan",
        "vol",
        "glyco_sites",
        "startAA",
        "endAA",
        "Seq_with_mod",
        "Glycopeptide_identifier",
        "Oxonium_ions",
        "bare_b_ions",
        "total_b_ions_possible",
        "bare_y_ions",
        "total_y_ions_possible",
        "b_ions_with_HexNAc",
        "y_ions_with_HexNAc",
        "b_ion_coverage",
        "y_ion_coverage",
        "Stub_ions",
        "percent_b_ion_with_HexNAc_coverage",
        "percent_y_ion_with_HexNAc_coverage",
        "percent_b_ion_coverage",
        "percent_y_ion_coverage",
        "#_of_stubs_found",
        'scan_id']

    f = open(output_file + ".csv", 'wb')
    dict_writer = csv.DictWriter(f, keys)
    dict_writer.writer.writerow(keys)
    dict_writer.writerows(scored_results)
    f.close()

    return output_file + '.csv'


if __name__ == '__main__':
    import sys
    main(sys.argv[1])
