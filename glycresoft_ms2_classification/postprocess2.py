import json
import sqlitedict
from .structure import glycans as glycan_lib
from .utils import try_get_outfile
from .prediction_tools import PredictionResults


def main(db_file, output_file=None):
    '''
    Scores observed fragment to theoretical ion pairings, eliminates ambiguity using tandem MS2 data
    :param db_file: Path to file created by match_ions.py, matching theoretical ions to observed spectra
    :param output_file: Path to file to write results to. Defaults to `db_file` + ".processed"
    '''
    if output_file is None:
        output_file = try_get_outfile(db_file, "processed")
    matched_ions = sqlitedict.SqliteDict(db_file, tablename="matched_ions")
    scored_results = []
    did_not_match = []
    metadata = sqlitedict.SqliteDict(db_file, tablename="metadata")
    metadata['db_file_name'] = db_file
    glycan_identities = metadata["glycan_identities"]

    for row in matched_ions.values():
        total_bare_b = float(row["total_b_ions"])

        total_bare_y = float(row["total_y_ions"])

        all_b_found = (row["b_ion_coverage"])
        num_b_found = float(len(all_b_found))
        percent_b = (num_b_found / total_bare_b)

        all_y_found = row["y_ion_coverage"]
        num_y_found = float(len(all_y_found))
        percent_y = (num_y_found / total_bare_y)

        b_HexNAc = row["b_ions_with_HexNAc"]
        b_HexNAc_found = len(b_HexNAc)
        total_possible_b_HexNAc = float(row["possible_b_ions_HexNAc"])
        if total_possible_b_HexNAc == 0:
            percent_b_HexNAc = 0
        else:
            percent_b_HexNAc = b_HexNAc_found / total_possible_b_HexNAc

        y_HexNAc = row["y_ions_with_HexNAc"]
        y_HexNAc_found = len(y_HexNAc)
        total_possible_y_HexNAc = float(
            row["possible_y_ions_HexNAc"])
        if total_possible_y_HexNAc == 0:
            percent_y_HexNAc = 0
        else:
            percent_y_HexNAc = y_HexNAc_found / total_possible_y_HexNAc

        stubs = row["Stub_ions"]

        Oxoniums = row["Oxonium_ions"]
        glycan_lib.parse_glycan_pack_string(row["Glycan"], glycan_identities)
        glycans = glycan_lib.parse_glycan_pack_string(row["Glycan"], glycan_identities)
        NeuAc = glycans.get("NeuAc", 0)

        bad_oxonium_ions = []

        if NeuAc != 0:
            pass

        else:
            try:
                bad_oxonium_ions.extend(d for d in Oxoniums if d['key'] in {"NeuAC", "NeuAc-H2O"})
            except:
                bad_oxonium_ions.extend(d for d in Oxoniums if d['key'] == 274.092 or d["key"] == 292.1026)
        if (all_y_found + all_b_found) != 0:
            scored_results.append({"MS1_Score": row["MS1_Score"],
                                   "Obs_Mass": row["Obs_Mass"],
                                   "Calc_mass": row["Calc_mass"],
                                   "ppm_error": row["ppm_error"],
                                   "Peptide": row["Peptide"],
                                   "Peptide_mod": row["Peptide_mod"],
                                   "Glycan": row["Glycan"],
                                   "vol": row["vol"],
                                   "glyco_sites": row["glyco_sites"],
                                   "startAA": row["startAA"],
                                   "endAA": row["endAA"],
                                   "Seq_with_mod": row["Seq_with_mod"],
                                   "Glycopeptide_identifier": row["Seq_with_mod"] + row["Glycan"],
                                   "Oxonium_ions": Oxoniums,
                                   "bare_b_ions": row["bare_b_ions"],
                                   "total_b_ions_possible": total_bare_b,
                                   "bare_y_ions": row["bare_y_ions"],
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
                                   "scan_id": row['scan_id'],
                                   "scan_id_range": row["scan_id_range"],
                                   "bad_oxonium_ions": bad_oxonium_ions,
                                   "glycan_composition": glycans
                                   })
        else:
            did_not_match.append(row["Glycopeptide_identifier"])

    results = {
        "metadata": dict(metadata),
        "predictions": scored_results,
        "intact_mass_no_fragments": did_not_match
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
