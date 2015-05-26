import functools
import csv
import itertools
import logging
import multiprocessing

import sqlitedict

import numpy as np

from glycresoft_ms2_classification.structure import glycans as glycan_lib
from glycresoft_ms2_classification.structure import modification
from glycresoft_ms2_classification.structure import composition

logger = logging.getLogger(__name__)

Modification = modification.Modification
ModificationTable = modification.ModificationTable
Composition = composition.Composition

Glycan = glycan_lib.Glycan


def glycosylate_callback(peptide, glycans, position_selector):
    n_sites = len(peptide.n_glycan_sequon_sites)
    result = []
    for glycan_count in range(1, n_sites + 1):
        glycan_cache = set()
        for glycans in itertools.product(glycans, repeat=glycan_count):
            # logger.debug("Mixing %s", glycans)
            # Check if we've already used this glycan combination before to avoid
            # duplicating effort for this peptide.
            if frozenset(glycans) in glycan_cache:
                continue
            else:
                glycan_cache.add(frozenset(glycans))
                # logger.debug("Cache Miss")
            # Consider all combinations of sites (order does not
            # matter)
            for sites in position_selector(peptide.n_glycan_sequon_sites, glycan_count):
                target = peptide.clone()
                glycan_composition = np.zeros(
                    len(glycans[0].composition))
                glycan_iter = iter(glycans)
                for site in sites:
                    glycan = glycan_iter.next()
                    hexnac = Modification("HexNAc")
                    hexnac.mass = glycan.mass
                    # There may already be a glycan included here by the proteomics
                    # search engine. Remove it  so we can do a general search.
                    if len(target[site][1]) > 0:
                        target.drop_modification(site, target[site][1][0])
                    target.add_modification(site, hexnac)
                    glycan_composition += np.array(
                        glycan.composition)
                glycan_composition_string = Glycan.glycan_composition_string(
                    glycan_composition)
                target.glycan = glycan_composition_string
                result.append(target)
    return result


class GlycopeptideHypothesis(object):
    def __init__(self, proteome=None, glycans=None):
        self.proteome = proteome
        self.peptides = []
        self.glycans = glycans

    def glycosylate_not_position_specific(self, peptide):
        n_sites = len(peptide.n_glycan_sequon_sites)
        for glycan_count in range(1, n_sites + 1):
            glycan_cache = set()
            for glycans in itertools.product(self.glycans, repeat=glycan_count):
                # Check if we've already used this glycan combination before to avoid
                # duplicating effort for this peptide.
                if frozenset(glycans) in glycan_cache:
                    continue
                else:
                    glycan_cache.add(frozenset(glycans))
                # Consider all combinations of sites (order does not matter)
                for sites in itertools.combinations(peptide.n_glycan_sequon_sites, glycan_count):
                    target = peptide.clone()
                    glycan_composition = np.zeros(len(glycans[0].composition))
                    glycan_iter = iter(glycans)
                    for site in sites:
                        glycan = glycan_iter.next()
                        hexnac = Modification("HexNAc")
                        hexnac.mass = glycan.mass
                        target.add_modification(site, hexnac)
                        glycan_composition += np.array(
                            glycan.composition)
                    glycan_composition_string = Glycan.glycan_composition_string(
                        glycan_composition)
                    target.glycan = glycan_composition_string
                    yield target

    def glycosylate_position_specific(self, peptide):
        n_sites = len(peptide.n_glycan_sequon_sites)
        for glycan_count in range(1, n_sites + 1):
            glycan_cache = set()
            for glycans in itertools.product(self.glycans, repeat=glycan_count):
                for sites in itertools.permutations(peptide.n_glycan_sequon_sites, glycan_count):
                    logger.debug(sites)
                    target = peptide.clone()
                    glycan_composition = np.zeros(
                        len(glycans[0].composition))
                    glycan_iter = iter(glycans)
                    for site in sites:
                        glycan = glycan_iter.next()
                        mod = glycan.as_modification()
                        target.add_modification(site, mod)
                        glycan_composition += np.array(
                            glycan.composition)
                    glycan_composition_string = "[{0}]".format(
                        ";".join(map(lambda x: str(int(x)), glycan_composition)))
                    target.glycan = glycan_composition_string
                    logger.debug(target)
                    yield target

    def combine(self, position_specific=False, include_unglycosylated_peptides=False):
        logger.debug("Combining glycans and peptides")
        if position_specific:
            position_selector = itertools.permutations
        else:
            position_selector = itertools.combinations
        for peptide in self.proteome.peptides(no_dups=True):
            logger.debug("Handling %s, %s glycosylation sites", peptide, peptide.n_glycan_sequon_sites)
            n_sites = len(peptide.n_glycan_sequon_sites)
            if n_sites == 0:
                if include_unglycosylated_peptides:
                    yield peptide.clone()
                continue
            if not position_specific:
                for glycan_count in range(1, n_sites + 1):
                    glycan_cache = set()
                    for glycans in itertools.product(self.glycans, repeat=glycan_count):
                        logger.debug("Mixing %s", glycans)
                        # Check if we've already used this glycan combination before to avoid
                        # duplicating effort for this peptide.
                        if frozenset(glycans) in glycan_cache:
                            continue
                        else:
                            glycan_cache.add(frozenset(glycans))
                            logger.debug("Cache Miss")
                        # Consider all combinations of sites (order does not
                        # matter)
                        for sites in position_selector(peptide.n_glycan_sequon_sites, glycan_count):
                            target = peptide.clone()
                            glycan_composition = np.zeros(
                                len(glycans[0].composition))
                            glycan_iter = iter(glycans)
                            for site in sites:
                                glycan = glycan_iter.next()
                                hexnac = Modification("HexNAc")
                                hexnac.mass = glycan.mass
                                # There may already be a glycan included here by the proteomics
                                # search engine. Remove it  so we can do a general search.
                                if len(target[site][1]) > 0:
                                    target.drop_modification(site, target[site][1][0])
                                target.add_modification(site, hexnac)
                                glycan_composition += np.array(
                                    glycan.composition)
                            glycan_composition_string = Glycan.glycan_composition_string(
                                glycan_composition)
                            target.glycan = glycan_composition_string
                            yield target
            else:
                for glycan_count in range(1, n_sites + 1):
                    glycan_cache = set()
                    for glycans in itertools.product(self.glycans, repeat=glycan_count):
                        for sites in position_selector(peptide.n_glycan_sequon_sites, glycan_count):
                            logger.debug(sites)
                            target = peptide.clone()
                            glycan_composition = np.zeros(
                                len(glycans[0].composition))
                            glycan_iter = iter(glycans)
                            for site in sites:
                                glycan = glycan_iter.next()
                                mod = glycan.as_modification()
                                target.add_modification(site, mod)
                                glycan_composition += np.array(
                                    glycan.composition)
                            glycan_composition_string = "[{0}]".format(
                                ";".join(map(lambda x: str(int(x)), glycan_composition)))
                            target.glycan = glycan_composition_string
                            logger.debug(target)
                            yield target

    def main(self, position_specific=False):
        for glycopeptide in self.combine():
            print(glycopeptide)
            prompt = raw_input()
            if prompt == "q":
                break

    def to_legacy(self, stream, protease=None, n_processes=2):
        legacy_format(
            stream, self.combine(), self.glycans[0].glycan_identities, protease)

    def pcombine(self, n_processes=2):
        position_selector = itertools.combinations
        task_fn = functools.partial(glycosylate_callback, glycans=self.glycans, position_selector=position_selector)
        pool = multiprocessing.Pool(n_processes, maxtasksperchild=1)
        for result in (pool.imap_unordered(task_fn, (p for p in self.proteome.peptides() if len(p.n_glycan_sequon_sites) > 0), chunksize=1)):
            for item in result:
                logger.debug("Finished %s", item)
                yield item
            del result


def legacy_format(stream, glycopeptides_iter, glycan_identities=None, protease=None):
    if glycan_identities is None:
        glycan_identities = []
    columns = ["Molecular Weight", "C", "Compositions"] + glycan_identities +\
              ["Adduct/Replacement", "Adduct Amount", "Peptide Sequence", "Peptide Modification",
               "Peptide Missed Cleavage Number", "Number of Glycan Attachment to Peptide", "Start AA",
               "End AA", "ProteinID"]

    handle = csv.writer(stream)
    handle.writerow(columns)
    for glycopeptide in glycopeptides_iter:
        line_data = map(str, [glycopeptide.mass, 0, glycopeptide.glycan] + glycopeptide.glycan[1:-1].split(";") +
                        ["/0", 0, str(glycopeptide), "/0",
                         0 if protease is None else protease.count(str(glycopeptide)),
                         glycopeptide.mod_index["HexNAc"], glycopeptide.start, glycopeptide.end,
                         glycopeptide.protein_id])

        handle.writerow(line_data)
