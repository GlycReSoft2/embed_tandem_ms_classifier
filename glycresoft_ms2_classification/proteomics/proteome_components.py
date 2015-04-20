import json

from collections import defaultdict
from itertools import ifilter


from ..structure import sequence
from ..utils import collectiontools

Sequence = sequence.Sequence


class Proteome(dict):
    def __init__(self, proteins_dict=None):
        if proteins_dict is None:
            proteins_dict = dict()
        dict.__init__(self, **proteins_dict)
        self.index = None

    def __iter__(self):
        for protein in self.values():
            yield protein

    def compress(self):
        for protein in self:
            protein.compress()

    def peptide_index(self):
        self.index = index_shared_peptides(self)
        return self.index

    def peptides(self, no_dups=True, protein_filter=lambda x: True, peptide_filter=lambda x: True):
        if no_dups:
            track = set()
            for peptide in self.peptides(no_dups=False, protein_filter=protein_filter, peptide_filter=peptide_filter):
                k = hash(peptide)
                if k not in track:
                    track.add(k)
                    yield peptide
        else:
            for protein in ifilter(protein_filter, self):
                for peptide in ifilter(peptide_filter, protein):
                    yield peptide


class ReferenceProtein(object):
    def __init__(self, sequence_dict):
        try:
            self.accession = sequence_dict["accession"]
            self.sequence = sequence_dict["Seq"]
            self.peptides = []
            self.metadata = {}
        except:
            print(sequence_dict)
            raise

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        rep = '''>{accession} with {num_peptides} peptides
{sequence}
'''.format(accession=self.accession, sequence=self.sequence,
           num_peptides=len(self.peptides))
        return rep

    def __getitem__(self, indices):
        return self.sequence[indices]

    def __hash__(self):
        return hash(self.accession)

    def sort_peptides(self):
        '''Group peptides by sequence features and sort them by peptide_score'''
        groups = collectiontools.groupby(self.peptides, lambda x: (str(x), (x.start, x.end), x.parent.accession))
        map(lambda x: x.sort(key=lambda y: y.peptide_score, reverse=True), groups.values())
        return groups

    def __iter__(self):
        '''Iterate over the best scoring instance of each unique peptide match'''
        groups = self.sort_peptides()
        peptides = map(lambda group: sorted(group, key=lambda x: x.peptide_score)[0], groups.values())
        for peptide in peptides:
            yield peptide

    @property
    def protein_score(self):
        for score_name in PROTEOMICS_SCORE:
            if score_name in self.metadata:
                return self.metadata[score_name]

    def compress(self):
        self.peptides = list(self)

    @property
    def n_glycan_sequon_sites(self):
        return sequence.find_n_glycosylation_sequons(self.sequence) +\
               self.metadata.get("n_glycan_sequon_sites", [])

    def describe(self):
        rep = '''>{accession} with {num_peptides} peptides
{sequence}
{metadata}
'''.format(accession=self.accession, sequence=self.sequence,
           num_peptides=len(self.peptides), metadata=json.dumps(self.metadata,
                                                                indent=2, sort_keys=True))
        return rep


def index_shared_peptides(reference_protein_list):
    index = defaultdict(list)
    for reference_protein in reference_protein_list:
        for peptide_match in reference_protein:
            index[peptide_match].append(
                (reference_protein.accession, peptide_match)
                )
    return index


PROTEOMICS_SCORE = ["PEAKS:peptideScore", "mascot:score", "PEAKS:proteinScore"]


class PeptideMatch(Sequence):

    def __init__(self, sequence, protein_id, start, end, is_decoy, **metadata):
        super(PeptideMatch, self).__init__(sequence)
        self.protein_id = protein_id
        self.start = start
        self.end = end
        self.is_decoy = is_decoy
        self.metadata = metadata
        self.parent = None
        self._base_sequence = sequence

    @property
    def n_glycan_sequon_sites(self):
        if self.parent is not None:
            return tuple(site - self.start for site in self.parent.n_glycan_sequon_sites
                         if self.start <= site < self.end and len(self.seq[site - self.start][1]) == 0)
        else:
            return sequence.find_n_glycosylation_sequons(self)

    def __hash__(self):
        return hash(str(self))

    def __repr__(self):
        seq = super(PeptideMatch, self).__str__()
        decoy = "(Decoy)" if self.is_decoy else ""
        rep = '''>{protein_id}[{start}:{end}] {decoy}
{seq}
{metadata}'''.format(protein_id=self.protein_id, start=self.start, end=self.end,
                     decoy=decoy, seq=seq, metadata=json.dumps(self.metadata, indent=2, sort_keys=True))
        return rep

    @property
    def peptide_score(self):
        for score_name in PROTEOMICS_SCORE:
            if score_name in self.metadata:
                return self.metadata[score_name]
