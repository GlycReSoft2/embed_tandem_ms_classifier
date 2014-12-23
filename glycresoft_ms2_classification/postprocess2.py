import json

from structure import glycans as glycan_lib
from utils import try_deserialize, try_get_outfile


def main(matched_ions_file, output_file=None):
    '''
    Scores observed fragment to theoretical ion pairings, eliminates ambiguity using tandem MS2 data
    :param matched_ions_file: Path to file created by match_ions.py, matching theoretical ions to observed spectra
    :param output_file: Path to file to write results to. Defaults to `matched_ions_file` + ".processed"
    '''
    if output_file is None:
        output_file = try_get_outfile(matched_ions_file, "processed")
    results = try_deserialize(matched_ions_file)
    scored_results = []
    metadata = results["metadata"]
    matched_ions = results["matched_ions"]

    glycan_identities = metadata["glycan_identities"]

    for rows in matched_ions:
        total_bare_b = float(rows["total_b_ions"])

        total_bare_y = float(rows["total_y_ions"])

        all_b_found = (rows["b_ion_coverage"])
        num_b_found = float(len(all_b_found))
        percent_b = (num_b_found / total_bare_b)

        all_y_found = rows["y_ion_coverage"]
        num_y_found = float(len(all_y_found))
        percent_y = (num_y_found / total_bare_y)

        b_HexNAc = rows["b_ions_with_HexNAc"]
        b_HexNAc_found = len(b_HexNAc)
        total_possible_b_HexNAc = float(rows["possible_b_ions_HexNAc"])
        if total_possible_b_HexNAc == 0:
            percent_b_HexNAc = 0
        else:
            percent_b_HexNAc = b_HexNAc_found / total_possible_b_HexNAc

        y_HexNAc = rows["y_ions_with_HexNAc"]
        y_HexNAc_found = len(y_HexNAc)
        total_possible_y_HexNAc = float(
            rows["possible_y_ions_HexNAc"])
        if total_possible_y_HexNAc == 0:
            percent_y_HexNAc = 0
        else:
            percent_y_HexNAc = y_HexNAc_found / total_possible_y_HexNAc

        stubs = rows["Stub_ions"]

        Oxoniums = rows["Oxonium_ions"]
        glycan_lib.parse_glycan_pack_string(rows["Glycan"], glycan_identities)
        glycans = glycan_lib.parse_glycan_pack_string(rows["Glycan"], glycan_identities)
        NeuAc = glycans.get("NeuAc", 0)

        bad_oxonium_ions = []

        if NeuAc != 0:
            pass

        else:
            bad_oxonium_ions.extend(
                d for d in Oxoniums if d['key'] == 274.092 or d["key"] == 292.1026)

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
                               "Oxonium_ions": Oxoniums,
                               "bare_b_ions": rows["bare_b_ions"],
                               "total_b_ions_possible": total_bare_b,
                               "bare_y_ions": rows["bare_y_ions"],
                               "total_y_ions_possible": total_bare_y,
                               "b_ions_with_HexNAc": b_HexNAc,
                               "y_ions_with_HexNAc": y_HexNAc,
                               "b_ion_coverage": all_b_found,
                               "y_ion_coverage": all_y_found,
                               "Stub_ions": stubs,
                               "percent_b_ion_with_HexNAc_coverage": percent_b_HexNAc,
                               "percent_y_ion_with_HexNAc_coverage": percent_y_HexNAc,
                               "percent_b_ion_coverage": percent_b,
                               "percent_y_ion_coverage": percent_y,
                               "scan_id": rows['scan_id'],
                               "scan_id_range": rows["scan_id_range"],
                               "bad_oxonium_ions": bad_oxonium_ions,
                               "glycan_composition": glycans
                               })

    results = {
        "metadata": metadata,
        "predictions": scored_results
    }
    if output_file:
        f = open(output_file + ".json", 'wb')
        json.dump(results, f)
        f.close()
        return (output_file + '.json', results)
    else:
        return None, results


def taskmain():
    import sys
    main(sys.argv[1])
