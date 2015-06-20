import re
import logging
from pyteomics import mzid
from ..structure import sequence
from ..structure import modification
from ..structure import residue

from .proteome_components import ReferenceProtein, PeptideMatch, Proteome

logger = logging.getLogger("mzid")

Sequence = sequence.Sequence
ModificationRule = modification.ModificationRule
Modification = modification.Modification
Residue = residue.Residue


MzIdentML = mzid.MzIdentML


class Parser(MzIdentML):
    def _retrieve_refs(self, info, **kwargs):
        """Retrieves and embeds the data for each attribute in `info` that
        ends in _ref. Removes the id attribute from `info`"""
        for k, v in dict(info).items():
            if k.endswith('_ref'):
                try:
                    info.update(self.get_by_id(v, retrieve_refs=True))
                    del info[k]
                    info.pop('id', None)
                except:
                    # logger.debug("%s not found", v)
                    info['skip'] = True


def extract_proteins(mzid_file):
    parser = Parser(mzid_file, retrieve_refs=True, iterative=False, build_id_cache=True)
    protein_iter = parser.iterfind("ProteinDetectionHypothesis", retrieve_refs=True, recursive=True, iterative=False)
    protein_index = Proteome({})
    # logger.debug("Begin iteration")
    for protein_dict in protein_iter:
        # logger.debug("Hit Protein")
        protein = handle_protein_hypothesis(protein_dict)
        protein_index[protein.accession] = protein
    return protein_index


def handle_protein_hypothesis(hypothesis_dict):
    protein = ReferenceProtein(hypothesis_dict)
    # logger.debug("Processing %s", protein.sequence)
    protein.metadata.update({k: v for k, v in hypothesis_dict.items()
                             if k not in {"PeptideHypothesis"}})
    protein.peptides = [peptide
                        for match in hypothesis_dict["PeptideHypothesis"]
                        for ident_item in match["SpectrumIdentificationItemRef"]
                        for peptide in convert_dict_to_sequence(ident_item, protein)
                        ]
    return protein


def convert_dict_to_sequence(sequence_dict, parent_protein):
    # logger.debug("Input: %r, Parent: %r", sequence_dict, parent_protein)
    base_sequence = sequence_dict["PeptideSequence"]
    peptide_sequence = Sequence(sequence_dict["PeptideSequence"])
    if "SubstitutionModification" in sequence_dict:
        subs = sequence_dict["SubstitutionModification"]
        for sub in subs:
            pos = sub['location'] - 1
            replace = Residue(sub["replacementResidue"])
            peptide_sequence[pos][0] = replace

    if "Modification" in sequence_dict:
        mods = sequence_dict["Modification"]
        for mod in mods:
            pos = mod["location"] - 1
            modification = Modification(mod["name"])
            if pos == -1:
                peptide_sequence.n_term = modification
            elif pos == len(peptide_sequence):
                peptide_sequence.c_term = modification
            else:
                peptide_sequence.add_modification(pos, modification)
    for evidence in sequence_dict["PeptideEvidenceRef"]:
        start = evidence["start"] - 1
        end = evidence["end"]
        found = parent_protein.sequence.find(base_sequence)
        if found == -1:
            raise ValueError("Peptide not found in Protein")
        if found != start:
            # logger.debug(
            #     "%r: %d <- %d, %s, %d", evidence["PeptideSequence"], found, start,
            #     parent_protein.accession, parent_protein.sequence.count(base_sequence))
            start = found
            end = start + len(base_sequence)
        is_decoy = evidence["isDecoy"]
        try:
            match = PeptideMatch(str(peptide_sequence), parent_protein.accession, start,
                                 end, is_decoy, **{k: v for k, v in sequence_dict.items() if k not in
                                                   exclude_keys_from_sequence_dict})
            match.parent = parent_protein
            # logger.debug("Produce: %r", match)
            yield match
        except:
            print(evidence)
            raise
exclude_keys_from_sequence_dict = set(("accession", "isDecoy", "PeptideEvidenceRef"))
