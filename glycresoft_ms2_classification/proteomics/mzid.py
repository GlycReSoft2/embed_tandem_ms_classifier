from pyteomics import mzid


from ..structure import sequence
from ..structure import modification

from .proteome_components import ReferenceProtein, PeptideMatch, Proteome

Sequence = sequence.Sequence
ModificationRule = modification.ModificationRule
Modification = modification.Modification


def extract_proteins(mzid_file):
    protein_iter = mzid.iterfind(mzid_file, "ProteinDetectionHypothesis", retrieve_refs=True, iterative=False)
    protein_index = Proteome({})
    for protein_dict in protein_iter:
        protein = handle_protein_hypothesis(protein_dict)
        protein_index[protein.accession] = protein
    return protein_index


def handle_protein_hypothesis(hypothesis_dict):
    protein = ReferenceProtein(hypothesis_dict)
    protein.metadata.update({k: v for k, v in hypothesis_dict.items()
                             if k not in {"PeptideHypothesis"}})
    protein.peptides = [peptide
                        for match in hypothesis_dict["PeptideHypothesis"]
                        for ident_item in match["SpectrumIdentificationItemRef"]
                        for peptide in convert_dict_to_sequence(ident_item, protein.accession)
                        ]
    return protein


def convert_dict_to_sequence(sequence_dict, protein_id=None):
    peptide_sequence = sequence_dict["PeptideSequence"]
    n_term = ""
    c_term = ""
    if "SubstitutionModification" in sequence_dict:
        sequence_list = list(peptide_sequence)
        subs = sequence_dict["SubstitutionModification"]
        for sub in subs:
            pos = sub['location'] - 1
            mod_str = "@{originalResidue}->{replacementResidue}".format(**sub)
            sequence_list[pos] = "{0}({1})".format(sub['replacementResidue'], mod_str)
        peptide_sequence = "".join(sequence_list)

    if "Modification" in sequence_dict:
        mods = sequence_dict["Modification"]
        sequence_list = list(peptide_sequence)
        for mod in mods:
            try:
                pos = mod["location"] - 1
                sequence_list[pos] += "({0})".format(mod["name"])
            except IndexError:
                if pos == -1:
                    n_term = "({0})-".format(mod["name"])
                elif pos == len(sequence_list):
                    c_term = "-({0})".format(mod["name"])
                else:
                    raise
        peptide_sequence = n_term + "".join(sequence_list) + c_term
    for evidence in sequence_dict["PeptideEvidenceRef"]:
        #evidence_id = evidence["id"]
        start = evidence["start"] - 1
        end = evidence["end"]
        is_decoy = evidence["isDecoy"]
        try:
            match = PeptideMatch(peptide_sequence, protein_id, start,
                                 end, is_decoy, **{k: v for k, v in sequence_dict.items() if k not in
                                                   exclude_keys_from_sequence_dict})
            yield match
        except:
            print(evidence)
            raise
exclude_keys_from_sequence_dict = set(("PeptideEvidenceRef", "accession", "start", "end",
                                       "isDecoy", "PeptideSequence", "Modification"))
