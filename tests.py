import os
import unittest
import warnings
import csv
import logging

logging.basicConfig(level=logging.DEBUG)
#multiprocessing_util.log_to_stderr()

from sqlitedict import SqliteDict

from glycresoft_ms2_classification.proteomics.msdigest_xml_parser import MSDigestParamters
from glycresoft_ms2_classification import theoretical_glycopeptide
from glycresoft_ms2_classification import calculate_fdr
from glycresoft_ms2_classification import entry_point
from glycresoft_ms2_classification import classify_matches

from glycresoft_ms2_classification.structure import sequence
from glycresoft_ms2_classification.structure import glycans
from glycresoft_ms2_classification.structure import modification
from glycresoft_ms2_classification.structure import composition

from glycresoft_ms2_classification.prediction_tools.false_discovery_rate import random_glycopeptide


class IonMatchingPipeline(unittest.TestCase):
    db_file_name = "test_data/USSR/Resultsof20131219_005.db"
    ms1_matching_output_file = "test_data/USSR/Resultsof20131219_005.csv"
    ms2_decon_file = "test_data/USSR/USSR_Grouped.yaml"
    glycosylation_sites_file = "test_data/USSR/USSR-glycosylation site list.txt"
    protein_prospector_file = "test_data/USSR/ProteinProspectorKK.xml"
    constant_modifications = ["Carbamidomethyl (C)"]
    variable_modifications = ["Deamidated (N)", "Deamidated (Q)"]
    method = "full_random_forest"
    methods = classify_matches.ModelTask.method_table.keys()
    ms1_match_tolerance = 1E-05
    ms2_match_tolerance = 2E-05
    num_procs = 6

    postprocessed_ions_file = "test_data/USSR/Resultsof20131219_005.processed.json"
    model_file_path = "test_data/USSR/Resultsof20131219_005.model.json"
    test_model_file_path = "test_data/USSRInfluenzaModel.json"
    classification_results_file = "test_data/USSR/Resultsof20131219_005.scored.json"
    model_eval_file = None

    def test_1_theoretical_ion_space_step(self):
        print("test_1_theoretical_ion_space_step")
        ms_digest = MSDigestParamters.parse(self.protein_prospector_file)
        theo_ions = entry_point.generate_theoretical_ion_space(
            self.ms1_matching_output_file, self.glycosylation_sites_file,
            ms_digest.constant_modifications, ms_digest.variable_modifications,
            ms_digest.enzyme, self.num_procs)
        self.assertTrue(os.path.exists(theo_ions))
        self.theoretical_ion_space_file = theo_ions
        theoretical_ions = SqliteDict(theo_ions, tablename="theoretical_search_space")
        sequence_set = theoretical_ions.itervalues()
        peptide_sequences = [
            sequence.Sequence(s["Seq_with_mod"]) for s in sequence_set]
        peptide_mods = set()
        for seq in peptide_sequences:
            for resid, mod in seq:
                peptide_mods.update((m.rule for m in mod))
        print(peptide_mods)

    def test_2_match_ions_step(self):
        print("test_2_match_ions_step")
        matches = entry_point.match_deconvoluted_ions(
            self.db_file_name, self.ms2_decon_file,
            self.ms1_match_tolerance, self.ms2_match_tolerance, self.num_procs)
        self.assertTrue(os.path.exists(matches))
        self.ms2_match_file = matches
        print(self.ms2_match_file)

    def test_3_postprocess_matches_step(self):
        print("test_3_postprocess_matches_step")
        self.postprocessed_ions_file = entry_point.postprocess_matches(
            self.db_file_name)
        self.assertTrue(os.path.exists(self.postprocessed_ions_file))
        print(self.postprocessed_ions_file)

    def test_4_build_model_step(self):
        print("test_4_build_model_step")
        self.model_file_path = entry_point.prepare_model_file(self.postprocessed_ions_file,
                                                              method=self.method)
        print(self.model_file_path)
        self.assertTrue(os.path.exists(self.model_file_path))

    def test_5_classify_with_model_step(self):
        print("test_5_classify_with_model_step")
        warnings.simplefilter(action="error")
        print(self.test_model_file_path)
        self.classification_results_file = entry_point.classify_data_by_model(self.postprocessed_ions_file,
                                                                              self.test_model_file_path,
                                                                              method=self.method)
        self.assertTrue(os.path.exists(self.classification_results_file))

    # def test_6_evaluate_model_step(self):
    #     print("test_6_evaluate_model_step")
    #     for method in self.methods:
    #         try:
    #             self.model_eval_file = entry_point.ModelDiagnosticsTask(
    #                 self.test_model_file_path, method).run()
    #             self.assertTrue(os.path.exists(self.model_eval_file))
    #         except IOError, e:
    #             # Windows doesn't like really long file names.
    #             print(e)

    def test_7_calculate_fdr_step(self):
        print("test_7_calculate_fdr_step")

        predicates = calculate_fdr.default_predicates()
        self.fdr_results = calculate_fdr.main(self.classification_results_file, self.ms2_decon_file,
                                              self.test_model_file_path, suffix_len=1,
                                              predicate_fns=predicates, n_processes=self.num_procs)
        self.assertTrue(os.path.exists(self.classification_results_file[:-5] + "_fdr.json"))

    # def test_8_apply_to_different_dataset(self):
    #     reclassified_file = entry_point.CompareModelsDiagnosticTask(self.test_model_file_path,
    #                                                                 self.extern_eval_file,
    #                                                                 method=self.method).run()


class TestTheoreticalIonSpaceProgram(unittest.TestCase):
    glycan_identities = [
        'GalNAcS', 'GalP', 'GalS', 'NeuGc', 'Hex', 'Pen', 'Fuc', 'Neu',
        'HexNAc', 'NeuAcAc', 'ManP', 'Kdn', 'HexN', 'HexA', 'GlcNAcS',
        'Xxx', 'GlcAS', 'GalNAcS2', 'Rha', 'Xyl', 'NeuAc', 'Water']

    def glycan_identity_extraction(self):
        return True  # The glycan identities need to be re-ordered
        result_file = "test_data/MS1-matching-output 20131219_005.csv"
        compo_dict = csv.DictReader(open(result_file, "r"), delimiter=",")

        colnames = compo_dict.fieldnames
        glycan_identity = theoretical_glycopeptide.get_glycan_identities(
            colnames)
        self.assertTrue(all([glycan_identity[i] == self.glycan_identities[i]
                             for i in range(len(self.glycan_identities))]))

sequence_str = "QQQHLFGSNVTDC(Carbamidomethyl)SGNFC(Carbamidomethyl)LFR"

b_ion_fragment_masses = [
    129.065853575,
    257.124431115,
    385.183008655,
    522.24192053,
    635.325984545,
    782.39439849,
    839.415862225,
    926.44789066,
    1040.49081813,
    1139.55923207,
    1240.60691058,
    1355.63385365,
    1515.66450215,
    1602.69653058,
    1659.71799432,
    1773.76092179,
    1920.82933573,
    2080.85998424,
    2193.94404825,
    2341.0124622,
]


class TestSequenceFragmentation(unittest.TestCase):

    def test_fragment_mass_calculations(self):
        seq_obj = sequence.Sequence(sequence_str)
        for i, frag in enumerate(seq_obj.get_fragments("B")):
            self.assertAlmostEqual(frag[0].mass, b_ion_fragment_masses[i])


class TestSequenceParsing(unittest.TestCase):
    seq_no_glycosites = "PEPTIDE"
    seq_1_glycosites = "PEPTINETIDE"
    seq_2_glycosites = "PEPTINNSTIDE"

    def test_n_glycan_sequon_finder(self):
        no_sites = sequence.find_n_glycosylation_sequons(self.seq_no_glycosites)
        self.assertTrue(len(no_sites) == 0)
        one_site = sequence.find_n_glycosylation_sequons(self.seq_1_glycosites)
        self.assertTrue(len(one_site) == 1)
        self.assertTrue(one_site[0] == 5)
        two_site = sequence.find_n_glycosylation_sequons(self.seq_2_glycosites)
        self.assertTrue(len(two_site) == 2)
        self.assertTrue(two_site[0] == 5)
        self.assertTrue(two_site[1] == 6)


class TestRandomSequenceGenerator(unittest.TestCase):
    terminals = ["K", "R"]
    glycans = glycans.load_from_file()[:50]
    const_modifications = [modification.Modification("Carbamidomethyl")]
    var_modifications = [modification.Modification("Deamidation"), modification.Modification("Oxidation")]
    tolerance = 10e-6

    def test_generate(self):
        target_mass = 4388.827053
        builder = random_glycopeptide.RandomGlycopeptideBuilder(tolerance=self.tolerance, )

class TestComposition(unittest.TestCase):

    def test_electron_mass(self):
        electron = composition.Composition("e")
        self.assertAlmostEqual(electron.mass, electron_mass)

    def test_proton_mass(self):
        proton = composition.Composition("p")
        self.assertAlmostEqual(proton.mass, proton_mass)

    def test_protein_mass(self):
        seq_obj = sequence.Sequence(sequence_str)
        self.assertAlmostEqual(seq_obj.mass, sequence_mass)

sequence_mass = 2514.11686191
electron_mass = 0.000549
proton_mass = 1.00727647

if __name__ == '__main__':
    unittest.main()
