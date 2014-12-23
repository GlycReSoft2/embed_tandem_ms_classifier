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
constants.PARTIAL_HEXNAC_LOSS = False
constants.EXCLUDE_B1 = True
