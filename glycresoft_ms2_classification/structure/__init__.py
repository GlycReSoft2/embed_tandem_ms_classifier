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
from argparse import Namespace

constants = Namespace()
## Constants
# Tokenizer Constants
constants.MOD_BEGIN = "("
constants.MOD_END = ")"
constants.GLYCAN_BEGIN = "["
constants.GLYCAN_END = "]"
# Sequence Fragment Constants
constants.FRAG_OFFSET = 1
constants.PARTIAL_HEXNAC_LOSS = True
constants.EXCLUDE_B1 = True


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
