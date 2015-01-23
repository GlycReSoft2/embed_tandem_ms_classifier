import itertools
from glycresoft_ms2_classification.ms import spectra

from collections import Counter, defaultdict

ObservedPrecursorSpectrum = spectra.ObservedPrecursorSpectrum
ObservedTandemSpectrum = spectra.ObservedTandemSpectrum


try:
    from cmath import fabs
except:
    from math import fabs

cdef double Proton = 1.007276035

cdef inline double ppm_error(double query, double match):
    return (match - query)/match

cpdef list match_observed_wrap(dict theoretical, str conn_string, double ms1_tol, double ms2_tol):
        db = spectra.MSMSSqlDB(conn_string)
        cdef list real_ions
        obs_mass = theoretical["Calc_mass"]
        real_ions = [ObservedPrecursorSpectrum.from_sql(r, db)
                     for r in db.ppm_match_tolerance_search(obs_mass, ms1_tol)]
        return match_observed_to_theoretical(theoretical, list(real_ions), ms2_tol)

cpdef list match_observed_to_theoretical(dict theoretical, list observed, double ms2_tolerance):
    cdef int scan_id, tandem_ms_ind

    annotate = defaultdict(list)
    did_match_counter = Counter()

    cdef list results
    results = []
    cdef double obs_mass, deprotonated_mass
    cdef int total_b, total_y, total_b_hex, total_y_hex
    cdef str theoretical_sequence
    theoretical_sequence = theoretical["Glycopeptide_identifier"]
    cdef list theoretical_ions, oxonoiums, b_type, all_b_ions
    cdef list b_HexNAc_type, y_type, all_y_ions, y_HexNAc_type
    cdef list stub_type

    cdef dict theo_ion


    cdef int i, tand_i, theo_i
    for i in range(len(observed)):
        real_ion = observed[i]
        scan_id = real_ion.scan_ids[0]
        tandem_ms_ind = real_ion.id
        did_match_counter[tandem_ms_ind] += 1

        oxoniums = list()
        theoretical_ions = theoretical["Oxonium_ions"]
        for theo_i in range(len(theoretical_ions)):
            theo_ion = theoretical_ions[theo_i]
            deprotonated_mass = theo_ion['mass'] - Proton
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    oxoniums.append({
                        "ion": obs_mass + Proton,
                        "ppm_error": ppm * 1e6,
                        "key": theo_ion['key'],
                        "scan_id": tandem_ms_ind
                    })
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion['key']))

        b_type = []
        all_b_ions = []

        theoretical_ions = theoretical["bare_b_ions"]
        total_b = len(theoretical_ions)

        for theo_i in range(len(theoretical_ions)):
            theo_ion = theoretical_ions[theo_i]
            deprotonated_mass = theo_ion['mass'] - Proton
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    b_type.append({
                        "obs_ion": obs_mass + Proton,
                        "key": theo_ion['key'],
                        "ppm_error": ppm * 1e6,
                        "scan_id": tandem_ms_ind
                    })
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion['key']))
        all_b_ions.extend(b_type)

        theoretical_ions = theoretical["b_ions_with_HexNAc"]
        total_b_hex = len(theoretical_ions)
        b_HexNAc_type = []
        for theo_i in range(len(theoretical_ions)):
            theo_ion = theoretical_ions[theo_i]
            deprotonated_mass = theo_ion['mass'] - Proton
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    b_HexNAc_type.append({
                        "obs_ion": obs_mass + Proton,
                        "key": theo_ion['key'],
                        "ppm_error": ppm * 1e6,
                        "scan_id": tandem_ms_ind
                    })
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion['key']))
        all_b_ions.extend(b_HexNAc_type)

        theoretical_ions = theoretical["bare_y_ions"]
        y_type = []
        all_y_ions = []
        total_y = len(theoretical_ions)
        for theo_i in range(len(theoretical_ions)):
            theo_ion = theoretical_ions[theo_i]
            deprotonated_mass = theo_ion['mass'] - Proton
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    y_type.append({
                        "obs_ion": obs_mass + Proton,
                        "key": theo_ion['key'],
                        "ppm_error": ppm * 1e6,
                        "scan_id": tandem_ms_ind
                    })
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion['key']))
        all_y_ions.extend(b_HexNAc_type)

        theoretical_ions = theoretical['y_ions_with_HexNAc']
        y_HexNAc_type = []
        total_y_hex = len(theoretical_ions)
        for theo_i in range(len(theoretical_ions)):
            theo_ion = theoretical_ions[theo_i]
            deprotonated_mass = theo_ion['mass'] - Proton
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    y_HexNAc_type.append({
                        "obs_ion": obs_mass + Proton,
                        "key": theo_ion['key'],
                        "ppm_error": ppm * 1e6,
                        "scan_id": tandem_ms_ind
                    })
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion['key']))
        all_y_ions.extend(b_HexNAc_type)

        theoretical_ions = theoretical['pep_stub_ions']
        stub_type = []
        for theo_i in range(len(theoretical_ions)):
            theo_ion = theoretical_ions[theo_i]
            deprotonated_mass = theo_ion['mass'] - Proton
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    stub_type.append({
                        "obs_ion": obs_mass + Proton,
                        "key": theo_ion['key'],
                        "ppm_error": ppm * 1e6,
                        "scan_id": tandem_ms_ind
                    })
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion['key']))



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
            "Oxonium_ions": oxoniums, "bare_b_ions": b_type, "possible_b_ions_HexNAc": total_b_hex,
            "total_b_ions": total_b, "bare_y_ions": y_type, "possible_y_ions_HexNAc": total_y_hex,
            "total_y_ions": total_y, "b_ions_with_HexNAc": b_HexNAc_type,
            "y_ions_with_HexNAc": y_HexNAc_type,
            "b_ion_coverage": all_b_ions, "y_ion_coverage": all_y_ions, "Stub_ions": stub_type,
            "scan_id": scan_id,
            "scan_id_range": []
        })
    return [results, did_match_counter, annotate]

