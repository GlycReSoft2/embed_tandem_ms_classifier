import os
import unittest

import csv

import theoretical_glycopeptide
import entry_point

from structure import sequence
from structure import composition


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

sequence_mass = 2514.11686191
electron_mass = 0.000549
proton_mass = 1.00727647


class TestIonMatchingPipelineProgram(unittest.TestCase):
    ms1_matching_output_file = "test_data/MS1-matching-output 20131219_005.csv"
    ms2_decon_file = "test_data/YAML-input-for-MS2-20131219_005.mzML.results"
    glycosylation_sites_file = "test_data/USSR-glycosylation site list.txt"
    constant_modifications = ["Carbamidomethyl (C)"]
    variable_modifications = ["Deamidated (N)", "Deamidated (Q)"]
    method = "full"
    ms1_match_tolerance = 1E-05
    ms2_match_tolerance = 2E-05

    theoretical_ion_space_file = "test_data/MS1-matching-output 20131219_005.theoretical_ions.csv"
    ms2_match_file = "test_data/MS1-matching-output 20131219_005.match_frags.csv"
    postprocessed_ions_file = "test_data/MS1-matching-output 20131219_005.processed.csv"
    model_file_path = "test_data/MS1-matching-output 20131219_005.model.csv"
    test_model_file_path = "test_data/Test Model.csv"
    classification_results_file = None
    model_eval_file = None

    def test_1_theoretical_ion_space_step(self):
        print("test_1_theoretical_ion_space_step")
        theo_ions = entry_point.generate_theoretical_ion_space(
            self.ms1_matching_output_file, self.glycosylation_sites_file,
            self.constant_modifications, self.variable_modifications)
        self.assertTrue(os.path.exists(theo_ions))
        self.theoretical_ion_space_file = theo_ions

    def test_2_match_ions_step(self):
        print("test_2_match_ions_step")
        matches = entry_point.match_deconvoluted_ions(
            self.theoretical_ion_space_file, self.ms2_decon_file,
            self.ms1_match_tolerance, self.ms2_match_tolerance)
        self.assertTrue(os.path.exists(matches))
        self.ms2_match_file = matches

    def test_3_postprocess_matches_step(self):
        print("test_3_postprocess_matches_step")
        self.postprocessed_ions_file = entry_point.postprocess_matches(
            self.ms2_match_file)
        self.assertTrue(os.path.exists(self.postprocessed_ions_file))

    def test_4_build_model_step(self):
        print("test_4_build_model_step")
        self.model_file_path = entry_point.prepare_model_file(self.postprocessed_ions_file,
                                                              method=self.method)
        self.assertTrue(os.path.exists(self.model_file_path))

    def test_5_classify_with_model_step(self):
        print("test_5_classify_with_model_step")
        self.classification_results_file = entry_point.classify_data_by_model(self.postprocessed_ions_file,
                                                                              self.test_model_file_path,
                                                                              method=self.method)
        self.assertTrue(os.path.exists(self.classification_results_file))

    def test_6_evaluate_model_step(self):
        print("test_6_evaluate_model_step")
        self.model_eval_file = entry_point.ModelDiagnosticsTask(self.model_file_path, self.method).run()
        self.assertTrue(os.path.exists(self.model_eval_file))


class TestTheoreticalIonSpaceProgram(unittest.TestCase):
    glycan_identities = [
        'GalNAcS', 'GalP', 'GalS', 'NeuGc', 'Hex', 'Pen', 'Fuc', 'Neu',
        'HexNAc', 'NeuAcAc', 'ManP', 'Kdn', 'HexN', 'HexA', 'GlcNAcS',
        'Xxx', 'GlcAS', 'GalNAcS2', 'Rha', 'Xyl', 'NeuAc', 'Water']

    def test_glycan_identity_extraction(self):
        #print("Reading %s" % result_file)
        result_file = "test_data/MS1-matching-output 20131219_005.csv"
        compo_dict = csv.DictReader(open(result_file, "r"), delimiter=",")

        colnames = compo_dict.fieldnames
        glycan_identity = theoretical_glycopeptide.get_glycan_identities(
            colnames)
        self.assertTrue(all([glycan_identity[i] == self.glycan_identities[i]
                             for i in range(len(self.glycan_identities))]))


class TestSequenceFragmentation(unittest.TestCase):

    def test_fragment_mass_calculations(self):
        seq_obj = sequence.Sequence(sequence_str)
        for i, frag in enumerate(seq_obj.get_fragments("B")):
            self.assertAlmostEqual(frag[0].mass, b_ion_fragment_masses[i])


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

if __name__ == '__main__':
    unittest.main()
