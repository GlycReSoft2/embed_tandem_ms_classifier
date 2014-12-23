import itertools
import logging
logging.basicConfig(level="DEBUG")

import numpy as np

from glycresoft_ms2_classification.structure import glycans as glycan_lib
from glycresoft_ms2_classification.structure import modification

Modification = modification.Modification
ModificationTable = modification.ModificationTable

Glycan = glycan_lib.Glycan


class GlycopeptideHypothesis(object):
    def __init__(self, peptides=None, glycans=None):
        if peptides is None:
            peptides = []
        if glycans is None:
            glycans = []
        self.peptides = peptides
        self.glycans = glycans

    def combine(self, position_specific=False):
        if position_specific:
            position_selector = itertools.permutations
        else:
            position_selector = itertools.combinations
        for peptide in self.peptides:
            n_sites = len(peptide.n_glycan_sequon_sites)
            if not position_specific:
                for glycan_count in range(1, n_sites + 1):
                    glycan_cache = set()
                    for glycans in itertools.product(self.glycans, repeat=glycan_count):
                        # Check if we've already used this glycan combination before to avoid
                        # duplicating effort for this peptide.
                        if not position_specific:
                            if frozenset(glycans) in glycan_cache:
                                continue
                            else:
                                glycan_cache.add(frozenset(glycans))
                        # Consider all combinations of sites (order does not matter)
                        for sites in position_selector(peptide.n_glycan_sequon_sites, glycan_count):
                            target = peptide.clone()
                            glycan_composition = np.zeros(len(glycans[0].composition))
                            glycan_iter = iter(glycans)
                            for site in sites:
                                glycan = glycan_iter.next()
                                hexnac = Modification("HexNAc")
                                hexnac.mass = glycan.mass
                                target.add_modification(site, hexnac)
                                glycan_composition += np.array(glycan.composition)
                            glycan_composition_string = Glycan.glycan_composition_string(glycan_composition)
                            target.glycan = glycan_composition_string
                            yield target
            else:
                for glycan_count in range(1, n_sites + 1):
                    glycan_cache = set()
                    for glycans in itertools.product(self.glycans, repeat=glycan_count):
                        for sites in position_selector(peptide.n_glycan_sequon_sites, glycan_count):
                            logging.debug(sites)
                            target = peptide.clone()
                            glycan_composition = np.zeros(len(glycans[0].composition))
                            glycan_iter = iter(glycans)
                            for site in sites:
                                glycan = glycan_iter.next()
                                mod = glycan.as_modification()
                                target.add_modification(site, mod)
                                glycan_composition += np.array(glycan.composition)
                            glycan_composition_string = "[{0}]".format(
                                ";".join(map(lambda x: str(int(x)), glycan_composition)))
                            target.glycan = glycan_composition_string
                            logging.debug(target)
                            yield target

    def main(self, position_specific=False):
        for glycopeptide in self.combine():
            print(glycopeptide)
            prompt = raw_input()
            if prompt == "q":
                break
