__all__ = [
    "sequence",
    "modification",
    "composition",
    "fragment",
    "glycan_masses",
    "residue",
    "stub_glycopeptides",
    "sequence_space"
]


# A few base types for doing type-based behavior changes
class PeptideSequenceBase(object):
    '''
    A base type for classes describing peptide sequences, with or without modifiations
    '''
    def get_sequence(self, start=0):
        return self.getSequence(start)

    def getSequence(self, start=0):
        pass


class ModificationBase(object):
    '''
    A base type for classes describing peptide sequence modifications
    '''
    pass


class ResidueBase(object):
    '''
    A base type for classes describing amino acid residues
    '''
