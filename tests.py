import unittest

import csv

import theoretical_glycopeptide

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


class TestTheoreticalIonSpaceProgram(unittest.TestCase):
    glycan_identities = ['GalNAcS', 'GalP', 'GalS', 'NeuGc', 'Hex', 'Pen', 'Fuc', 'Neu',
                         'HexNAc', 'NeuAcAc', 'ManP', 'Kdn', 'HexN', 'HexA', 'GlcNAcS',
                         'Xxx', 'GlcAS', 'GalNAcS2', 'Rha', 'Xyl', 'NeuAc', 'Water']

    def test_glycan_identity_extraction(self):
        #print("Reading %s" % result_file)
        result_file = "test_data/ResultOfData_input_for_MS1_matching_20131219_005_isos.csv"
        compo_dict = csv.DictReader(open(result_file, "r"), delimiter=",")

        colnames = compo_dict.fieldnames
        glycan_identity = theoretical_glycopeptide.get_glycan_identities(colnames)
        self.assertTrue(all([glycan_identity[i] == self.glycan_identities[i]
                             for i in range(len(self.glycan_identities))]))


class TestSequenceFragmentation(unittest.TestCase):

    def test_fragment_mass_calculations(self):
        seq_obj = sequence.Sequence(sequence_str)
        for i, frag in enumerate(seq_obj.getFragments("B")):
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
