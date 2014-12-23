import json

from pyteomics import mzid

from ..structure import sequence
from ..structure import modification

Sequence = sequence.Sequence
ModificationRule = modification.ModificationRule
Modification = modification.Modification


class ReferenceProtein(object):
    def __init__(self, sequence_dict):
        try:
            self.accession = sequence_dict["accession"]
            self.sequence = sequence_dict["Seq"]
            self.is_decoy = sequence_dict["isDecoy"]
            self.peptides = []
        except:
            print(sequence_dict)
            raise

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        decoy = "" if not self.is_decoy else "(Decoy)"
        rep = '''>{accession} with {num_peptides} peptides {decoy}
{sequence}'''.format(accession=self.accession, sequence=self.sequence,
                     num_peptides=len(self.peptides), decoy=decoy)
        return rep

    def __hash__(self):
        return hash()


class MzIdentML(object):
    def __init__(self, file_path, record=True):
        self.file_path = file_path
        self.reader = mzid.read(open(self.file_path), recursive=True, retrieve_refs=True)
        self.proteins = dict()
        self.record = record

    def __iter__(self):
        generator = iter(self.reader)
        for element in generator:
            yield PeptideSpectrumMatchSet(self, self.record, **element)

    def main(self, include_decoys=False):
        for matches in self:
            for peptide in matches.peptide_matches():
                if include_decoys and peptide.is_decoy:
                    continue
                print("Number of reference proteins: {0}".format(len(self.proteins)))
                print(repr(peptide))
                p = raw_input()
                if p == "q":
                    return


def convert_dict_to_sequence(sequence_dict):
    peptide_sequence = sequence_dict.pop("PeptideSequence")
    if "Modification" in sequence_dict:
        mods = sequence_dict.get("Modification")
        sequence_list = list(peptide_sequence)
        for mod in mods:
            sequence_list[mod["location"] - 1] += "({0})".format(mod["name"])
        peptide_sequence = "".join(sequence_list)
    metadata = sequence_dict.pop("PeptideEvidenceRef")[0]
    accession = metadata.pop("accession")
    start = metadata.pop("start")
    end = metadata.pop("end")
    is_decoy = metadata.pop("isDecoy")
    return PeptideMatch(peptide_sequence, accession, start,
                        end, is_decoy, **sequence_dict)


class PeptideSpectrumMatchSet(dict):

    def __init__(self, source, record=True, **kwargs):
        super(PeptideSpectrumMatchSet, self).__init__(**kwargs)
        self.source = source
        self.record = record

    def peptide_matches(self):
        for peptide_match in self["SpectrumIdentificationItem"]:
            match_parent = peptide_match['PeptideEvidenceRef'][0]
            if match_parent["accession"] not in self.source.proteins:
                parent = ReferenceProtein(match_parent)
                self.source.proteins[match_parent["accession"]] = parent
            else:
                parent = self.source.proteins[match_parent["accession"]]
            peptide = convert_dict_to_sequence(peptide_match)
            if self.record:
                self.source.proteins[parent.accession].peptides.append(peptide)
            yield peptide


class PeptideMatch(Sequence):
    def __init__(self, sequence, parent_accession, start, end, is_decoy, **metadata):
        super(PeptideMatch, self).__init__(sequence)
        self.accession = parent_accession
        self.start = start
        self.end = end
        self.is_decoy = is_decoy
        self.metadata = metadata

    def __repr__(self):
        seq = super(PeptideMatch, self).__str__()
        decoy = "(Decoy)" if self.is_decoy else ""
        rep = '''>{accession}[{start}:{end}] {decoy}
{seq}
{metadata}'''.format(accession=self.accession, start=self.start, end=self.end,
                     decoy=decoy, seq=seq, metadata=json.dumps(self.metadata, indent=2, sort_keys=True))
        return rep
