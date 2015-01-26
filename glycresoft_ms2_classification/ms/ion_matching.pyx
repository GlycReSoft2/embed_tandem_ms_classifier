# cython: profile=True
cimport cython

from collections import Counter, defaultdict

from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

try:
    import ujson as json
except:
    import json

from libc.math cimport fabs
from glycresoft_ms2_classification.ms import spectra

cdef class TheoreticalIon:
    cdef double mass
    cdef str key

    def __init__(self, double mass, str key):
        self.mass = mass
        self.key = key

    def __repr__(TheoreticalIon self):
        return "<TheoreticalIon {0} {1}>".format(self.key, self.mass)


cdef list theoretical_ion_list(list dict_ions):
    cdef:
        list results
        int i
        double mass
        str key
        dict ion_d
        TheoreticalIon ion_s
    results = [0] * len(dict_ions)
    for i in range(len(dict_ions)):
        ion_d = dict_ions[i]
        mass = ion_d['mass']
        key = str(ion_d["key"])
        ion_s = TheoreticalIon(mass, key)
        results[i] = ion_s
    return results


cdef class CObservedPrecursorSpectrum:
    cdef public double neutral_mass
    cdef int charge
    cdef list scans
    cdef public list scan_ids
    cdef dict other_data
    cdef public list tandem_data
    cdef public int _iterkey

    def __init__(self, list scan_data, list scan_ids, int z, double neutral_mass, list tandem_data, **other_data):
        self.scans = scan_data
        self.scan_ids = scan_ids
        self.charge = z
        self.neutral_mass = neutral_mass
        self.tandem_data = tandem_data
        self.other_data = other_data
        self._iterkey = -1

    property id:
        def  __get__(self):
            return self._iterkey

        def __set__(self, int value):
            self._iterkey = value

@cython.boundscheck(False)
@cython.wraparound(False)
cdef CObservedPrecursorSpectrum CObservedPrecursorSpectrum_from_sql(object row):
    cdef:
        double neutral_mass
        list scans, scan_ids, tandem_data
        int charge, _iterkey
        dict other_datai
        CObservedPrecursorSpectrum inst
    scans = json.loads(row["scan_data"])
    scan_ids = [s['id'] for s in scans]
    neutral_mass = row['neutral_mass']
    charge = row['charge']
    tandem_data = [CObservedTandemSpectrum(jdict["neutral_mass"], jdict["charge"], jdict["intensity"],
                                           jdict["id"], jdict["annotation"])
                   for jdict in json.loads(row['tandem_data'])]
    other_data = json.loads(row['other_data'])

    inst = CObservedPrecursorSpectrum(scans, scan_ids, charge, neutral_mass, tandem_data, **other_data)
    inst._iterkey = row['precursor_id']

    return inst


cdef class CObservedTandemSpectrum:
    cdef public double mass
    cdef int charge
    cdef double intensity
    cdef public  int id
    cdef list annotation
    cdef dict other_data

    def __init__(self, double neutral_mass, int charge, double intensity, int id, list annotation, **data):
        self.mass = neutral_mass
        self.charge = charge
        self.intensity = intensity
        self.id = id
        if annotation is None:
            annotation = []
        self.annotation = annotation
        self.other_data = data

    property neutral_mass:
        def __get__(self):
            return self.mass

cdef double Proton = 1.007276035

cdef inline double ppm_error(double query, double match):
    return (query - match)/match

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef list match_observed_wrap(dict theoretical, str observed_ions_conn_string,
                               double ms1_tolerance, double ms2_tolerance):
        db = spectra.MSMSSqlDB(observed_ions_conn_string)
        cdef list real_ions
        obs_mass = theoretical["Calc_mass"]
        real_ions = [CObservedPrecursorSpectrum_from_sql(r)
                     for r in db.ppm_match_tolerance_search(obs_mass, ms1_tolerance)]
        return match_observed_to_theoretical(theoretical, list(real_ions), ms2_tolerance)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef list match_observed_to_theoretical(dict theoretical, list observed, double ms2_tolerance):
    cdef int scan_id, tandem_ms_ind

    annotate = defaultdict(list)
    did_match_counter = Counter()

    cdef list results
    results = []
    cdef double obs_mass, deprotonated_mass, ppm
    cdef int total_b, total_y, total_b_hex, total_y_hex
    cdef str theoretical_sequence
    theoretical_sequence = str(theoretical["Glycopeptide_identifier"])
    cdef list theoretical_ions, oxonoiums, b_type, all_b_ions
    cdef list b_HexNAc_type, y_type, all_y_ions, y_HexNAc_type
    cdef list stub_type

    cdef:
        list theo_oxoniums, theo_b_type, theo_b_HexNAc_type
        list theo_y_type, theo_y_HexNAc_type, theo_stubs
        TheoreticalIon theo_ion
        CObservedPrecursorSpectrum real_ion
        CObservedTandemSpectrum real_tandem

    theo_oxoniums = theoretical_ion_list(theoretical["Oxonium_ions"])
    theo_b_type = theoretical_ion_list(theoretical["bare_b_ions"])
    theo_b_HexNAc_type = theoretical_ion_list(theoretical["b_ions_with_HexNAc"])
    theo_y_type = theoretical_ion_list(theoretical["bare_y_ions"])
    theo_y_HexNAc_type = theoretical_ion_list(theoretical["y_ions_with_HexNAc"])
    theo_stubs = theoretical_ion_list(theoretical["pep_stub_ions"])

    cdef int i, tand_i, theo_i
    for i in range(len(observed)):
        real_ion = observed[i]
        scan_id = real_ion.scan_ids[0]
        tandem_ms_ind = real_ion.id
        did_match_counter[tandem_ms_ind] += 1

        oxoniums = list()
        for theo_i in range(len(theo_oxoniums)):
            theo_ion = theo_oxoniums[theo_i]
            deprotonated_mass = theo_ion.mass
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass + Proton, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    oxoniums.append({
                        "ion": obs_mass + Proton,
                        "ppm_error": ppm * 1e6,
                        "key": theo_ion.key,
                        "scan_id": tandem_ms_ind
                    })
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion.key))

        b_type = []
        all_b_ions = []

        total_b = len(theo_b_type)
        b_type = []
        all_b_ions = []

        for theo_i in range(len(theo_b_type)):
            theo_ion = theo_b_type[theo_i]
            deprotonated_mass = theo_ion.mass - Proton
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    b_type.append({
                        "obs_ion": obs_mass + Proton,
                        "key": theo_ion.key,
                        "ppm_error": ppm * 1e6,
                        "scan_id": tandem_ms_ind
                    })
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion.key))
        all_b_ions.extend(b_type)

        total_b_hex = len(theo_b_HexNAc_type)
        b_HexNAc_type = []
        for theo_i in range(len(theo_b_HexNAc_type)):
            theo_ion = theo_b_HexNAc_type[theo_i]
            deprotonated_mass = theo_ion.mass - Proton
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    b_HexNAc_type.append({
                        "obs_ion": obs_mass + Proton,
                        "key": theo_ion.key,
                        "ppm_error": ppm * 1e6,
                        "scan_id": tandem_ms_ind
                    })
                    all_b_ions.append(
                    {
                        "obs_ion": obs_mass + Proton,
                        "key": theo_ion.key.split("+")[0],
                        "ppm_error": ppm * 1e6,
                        "scan_id": tandem_ms_ind
                    }
                    )
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion.key))

        y_type = []
        all_y_ions = []
        total_y = len(theo_y_type)
        for theo_i in range(len(theo_y_type)):
            theo_ion = theo_y_type[theo_i]
            deprotonated_mass = theo_ion.mass - Proton
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    y_type.append({
                        "obs_ion": obs_mass + Proton,
                        "key": theo_ion.key,
                        "ppm_error": ppm * 1e6,
                        "scan_id": tandem_ms_ind
                    })
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion.key))

        total_y_hex = len(theo_y_HexNAc_type)
        y_HexNAc_type = []
        for theo_i in range(len(theo_y_HexNAc_type)):
            theo_ion = theo_y_HexNAc_type[theo_i]
            deprotonated_mass = theo_ion.mass - Proton
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    y_HexNAc_type.append({
                        "obs_ion": obs_mass + Proton,
                        "key": theo_ion.key,
                        "ppm_error": ppm * 1e6,
                        "scan_id": tandem_ms_ind
                    })
                    all_y_ions.append(
                    {
                        "obs_ion": obs_mass + Proton,
                        "key": theo_ion.key.split("+")[0],
                        "ppm_error": ppm * 1e6,
                        "scan_id": tandem_ms_ind
                    }
                    )
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion.key))

        stub_type = []
        for theo_i in range(len(theo_stubs)):
            theo_ion = theo_stubs[theo_i]
            deprotonated_mass = theo_ion.mass - Proton
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    stub_type.append({
                        "obs_ion": obs_mass + Proton,
                        "key": theo_ion.key,
                        "ppm_error": ppm * 1e6,
                        "scan_id": tandem_ms_ind
                    })
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion.key))

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

