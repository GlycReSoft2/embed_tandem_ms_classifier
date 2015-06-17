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
    '''
    Represents a single predicted ion
    '''
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


cpdef list merge_observed_matches(list matches):
    merge_dict = defaultdict(list)
    cdef int i = 0
    cdef Match match
    cdef list results = []
    for i in range(0, len(matches)):
        match = matches[i]
        merge_dict[match.key].append(match)

    for key, group in merge_dict.items():
        results.append(merge_group(group).as_dict())

    return results


cpdef MatchGroup merge_group(list matches):
    cdef:
        double best_ppm = fabs(matches[0].ppm_error)
        double current_ppm
        Match best_match = matches[0]
        Match match
        int i = 0
        double sum_intensity = 0
        dict scan_id_to_ppm = dict()
        MatchGroup merged

    scan_id_to_ppm[matches[0].scan_id] = matches[0].ppm_error

    for i in range(i, len(matches)):
        match = matches[i]
        current_ppm = fabs(match.ppm_error)
        if current_ppm < best_ppm:
            best_ppm = current_ppm
            best_match = match
        scan_id_to_ppm[match.scan_id] = dict(ppm_error=match.ppm_error, 
                                             intensity=match.intensity,
                                             charge=match.charge)
        sum_intensity += match.intensity

    merged = MatchGroup(best_match.ppm_error, best_match.obs_ion, sum_intensity,
                        best_match.key, best_match.scan_id, scan_id_to_ppm)
    return merged


cdef class MatchGroup:
    cdef:
        public double ppm_error
        public double obs_ion
        public str key
        public int scan_id
        public double intensity
        public dict scan_map

    def __init__(self, double ppm_error, double obs_ion, double intensity,
                 str key, int scan_id, dict scan_map):
        if scan_map is None:
            scan_map = dict()
        self.ppm_error = ppm_error
        self.obs_ion = obs_ion
        self.intensity = intensity
        self.key = key
        self.scan_id = scan_id
        self.scan_map = scan_map

    cpdef dict as_dict(self):
        cdef dict result = {}
        result['ppm_error'] = self.ppm_error
        result['obs_ion'] = self.obs_ion
        result["intensity"] = self.intensity
        result['key'] = self.key
        result['scan_id'] = self.scan_id
        result['scan_map'] = dict(self.scan_map)
        return result

    def __repr__(MatchGroup self):
        return "<MatchGroup {0} {1} {2} {3} {4}>".format(self.key, self.ppm_error, 
                                                         self.intensity, self.scan_id, self.scan_map)


cdef class Match:
    cdef:
        public double ppm_error
        public double obs_ion
        public double intensity
        public int charge
        public str key
        public int scan_id

    def __init__(self, double ppm_error, double obs_ion, double intensity,
                 int charge, str key, int scan_id=-1):
        self.ppm_error = ppm_error
        self.obs_ion = obs_ion
        self.intensity = intensity
        self.charge = charge
        self.key = key + ':' + str(charge)
        self.scan_id = scan_id

    cpdef dict as_dict(self):
        cdef dict result = {}
        result['ppm_error'] = self.ppm_error
        result['obs_ion'] = self.obs_ion
        result['intensity'] = self.intensity
        result['charge'] = self.charge
        result['key'] = self.key
        result['scan_id'] = self.scan_id
        return result

    def __repr__(Match self):
        return "<Match {0} {1} {2} {3} {4}>".format(self.key, self.ppm_error, self.intensity, 
                                                    self.charge, self.scan_id)

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
        dict other_data
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
    cdef public int charge
    cdef public double intensity
    cdef public int id
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
    cdef str theoretical_sequence
    theoretical_sequence = str(theoretical["Glycopeptide_identifier"])
    
    cdef int total_b = 0, total_y = 0, total_b_hex = 0, total_y_hex = 0
    cdef list theoretical_ions = [], oxonoiums = [], b_type = [], all_b_ions = []
    cdef list b_HexNAc_type = [], y_type = [], all_y_ions = [], y_HexNAc_type = []
    cdef list stub_type = [], scan_id_range = [], oxoniums = []

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
        scan_id_range.append(scan_id)
        tandem_ms_ind = real_ion.id
        did_match_counter[tandem_ms_ind] += 1


        for theo_i in range(len(theo_oxoniums)):
            theo_ion = theo_oxoniums[theo_i]
            deprotonated_mass = theo_ion.mass
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass + Proton, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    oxoniums.append(Match(
                        ppm * 1e6,
                        obs_mass + Proton,
                        real_tandem.intensity,
                        real_tandem.charge,
                        theo_ion.key,
                        real_tandem.id))
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion.key))


        total_b = len(theo_b_type)
        # b_type = []
        # all_b_ions = []

        for theo_i in range(len(theo_b_type)):
            theo_ion = theo_b_type[theo_i]
            deprotonated_mass = theo_ion.mass - Proton
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    b_type.append(Match(
                        ppm * 1e6,
                        obs_mass + Proton,
                        real_tandem.intensity,
                        real_tandem.charge,
                        theo_ion.key,
                        real_tandem.id
                    ))
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion.key))
        all_b_ions.extend(b_type)


        total_b_hex = len(theo_b_HexNAc_type)
        # b_HexNAc_type = []
        for theo_i in range(len(theo_b_HexNAc_type)):
            theo_ion = theo_b_HexNAc_type[theo_i]
            deprotonated_mass = theo_ion.mass - Proton
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    b_HexNAc_type.append(Match(
                        ppm * 1e6,
                        obs_mass + Proton,
                        real_tandem.intensity,
                        real_tandem.charge,
                        theo_ion.key,
                        real_tandem.id
                    ))
                    all_b_ions.append(Match(
                        ppm * 1e6,
                        obs_mass + Proton,
                        real_tandem.intensity,
                        real_tandem.charge,
                        theo_ion.key.split("+")[0],
                        real_tandem.id
                    ))
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion.key))


        total_y = len(theo_y_type)
        for theo_i in range(len(theo_y_type)):
            theo_ion = theo_y_type[theo_i]
            deprotonated_mass = theo_ion.mass - Proton
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    y_type.append(Match(
                        ppm * 1e6,
                        obs_mass + Proton,
                        real_tandem.intensity,
                        real_tandem.charge,
                        theo_ion.key,
                        real_tandem.id
                    ))
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion.key))


        total_y_hex = len(theo_y_HexNAc_type)

        for theo_i in range(len(theo_y_HexNAc_type)):
            theo_ion = theo_y_HexNAc_type[theo_i]
            deprotonated_mass = theo_ion.mass - Proton
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    y_HexNAc_type.append(Match(
                        ppm * 1e6,
                        obs_mass + Proton,
                        real_tandem.intensity,
                        real_tandem.charge,
                        theo_ion.key,
                        real_tandem.id
                    ))
                    all_y_ions.append(Match(
                        ppm * 1e6,
                        obs_mass + Proton,
                        real_tandem.intensity,
                        real_tandem.charge,
                        theo_ion.key.split("+")[0],
                        real_tandem.id
                    ))
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion.key))



        for theo_i in range(len(theo_stubs)):
            theo_ion = theo_stubs[theo_i]
            deprotonated_mass = theo_ion.mass - Proton
            for tand_i in range(len(real_ion.tandem_data)):
                real_tandem = real_ion.tandem_data[tand_i]
                obs_mass = real_tandem.mass
                ppm = ppm_error(obs_mass, deprotonated_mass)
                if fabs(ppm) <= ms2_tolerance:
                    stub_type.append(Match(
                        ppm * 1e6,
                        obs_mass + Proton,
                        real_tandem.intensity,
                        real_tandem.charge,
                        theo_ion.key,
                        real_tandem.id
                    ))
                    annotate[real_tandem.id].append((theoretical_sequence, theo_ion.key))

    if len(did_match_counter) > 0:
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
            
            "total_b_ions": total_b,
            "possible_b_ions_HexNAc": total_b_hex,
            "total_y_ions": total_y,
            "possible_y_ions_HexNAc": total_y_hex,

            "Oxonium_ions": merge_observed_matches(oxoniums),
            "bare_b_ions": merge_observed_matches(b_type),
            "bare_y_ions": merge_observed_matches(y_type),
            "b_ions_with_HexNAc": merge_observed_matches(b_HexNAc_type),
            "y_ions_with_HexNAc": merge_observed_matches(y_HexNAc_type),
            "b_ion_coverage": merge_observed_matches(all_b_ions), 
            "y_ion_coverage": merge_observed_matches(all_y_ions),
            "Stub_ions": merge_observed_matches(stub_type),
            
            "scan_id": scan_id_range[0],
            "scan_id_range": scan_id_range,
            "protein_id": theoretical.get("protein_id")
        })

    return [results, did_match_counter, annotate]

