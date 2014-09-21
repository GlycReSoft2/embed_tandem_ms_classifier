import copy
import warnings

from structure import PeptideSequenceBase
from structure.composition import Composition
from structure.fragment import Fragment
from structure.modification import Modification
from structure.residue import Residue


# sequence must be a string and not unicode
def sequence_tokenizer(sequence):
    state = "aa"
    chunks = []
    current_aa = ""
    current_mod = ""
    paren_level = 0
    i = 0
    while i < len(sequence):
        next_char = sequence[i]
        if next_char is "(":
            if state == "aa":
                state = "mod"
                assert paren_level == 0
                paren_level += 1
            elif state == "mod":
                paren_level += 1
                current_mod += next_char
        elif next_char == ")":
            if state == "aa":
                raise Exception("Invalid Sequence. ) found outside of modification.")
            elif state == "mod":
                paren_level -= 1
                if paren_level == 0:
                    state = 'aa'
                    new_chunk = [current_aa, current_mod]
                    chunks.append(new_chunk)
                    current_mod = ""
                    current_aa = ""
                else:
                    current_mod += next_char
        elif state == "aa":
            if(current_aa != ""):
                new_chunk = [current_aa, current_mod]
                chunks.append(new_chunk)
                current_mod = ""
                current_aa = ""
            current_aa += next_char
        elif state == "mod":
            current_mod += next_char
        else:
            raise Exception("Unknown Tokenizer State", current_aa, current_mod, i, next_char)
        i += 1
    if current_aa != "":
        new_chunk = [current_aa, ""]
        chunks.append(new_chunk)
    return chunks


def strip_modifications(sequence):
    tokens = sequence_tokenizer(sequence)
    return ''.join([residue for residue, mod in tokens])


def sequence_to_mass(sequence):
    '''Fast and unapologetic or validating sequence to mass'''
    mass = 0.0
    tokens = sequence_tokenizer(sequence)
    for residue, mod in tokens:
        mass += Residue.mass_by_name(residue)
        mass += Modification.mass_by_name(mod)
    return mass


class Sequence(PeptideSequenceBase):
    """Name and mass of residues (glycan/peptide)"""

    def __init__(self, sequence):
        seq_list = sequence_tokenizer(sequence)
        self.mass = 0.0
        self.seq = []
        for item in seq_list:
            res = Residue()
            res.by_symbol(item[0])
            if item[1] != '':
                try:
                    mod = Modification(item[1])
                except:
                    print(item, "Modification failed!", sequence, seq_list)
                    raise
                self.seq.append([res, [mod]])
                self.mass += mod.mass
            else:
                self.seq.append([res, []])
            self.mass += res.mass
        self.length = len(self.seq)

        self.mass += Composition("H2O").mass

    def __repr__(self):
        rep = "{seq}[{mass}]".format(**self.__dict__)
        return rep

    def __len__(self):
        return self.length

    def __iter__(self):
        '''
            Makes the sequence iterable
        '''
        for i in self.seq:
            yield(i)

    def __getitem__(self, index):
        '''
            A Pythonic way to access an index in the sequence
        '''
        return self.seq[index]

    def __setitem__(self, index, value):
        '''
            A Pythonic way to set the value at an index, but does not
            try to validate the result.
        '''
        self.seq[index] = value

    def add_mass(self, pos, mass):
        """
            This function will direcly modify the mass of the residue instead of appending the modificaion.
            Use it with caution!
        """
        self.seq[pos][0].mass += mass
        self.mass += mass

    def get_fragments(self, type):
        """ Return a list of mass values"""

        mass_shift = 0.0

        # The set of modification names.
        mod_dict = {}

        # The key is the position, the value is an array of fragments.
        # And the first element is always bare fragment.
        # The total number of HexNAc on the fragment should be recorded.\

        if type == 'B':
            seq_list = self.seq

            mass_shift = Composition('H').mass - Composition('e').mass

        elif type == 'Y':
            mass_shift = Composition('H2O').mass + Composition(
                'H').mass - Composition('e').mass
            seq_list = list(reversed(self.seq))

        # Updated!!! The modified version should include the
        # first residue and the complete sequence but with
        # varied number of modifications.
        # Bare mass.
        current_mass = 0
        hn_num = 0
        for idx in range(len(seq_list) - 1):
            for mod in seq_list[idx][1]:
                # HexNAc is treated differently.
                if mod.name == 'HexNAc':
                    hn_num += 1
                    continue
                if mod.name in mod_dict:
                    mod_dict[mod.name] += 1
                else:
                    mod_dict[mod.name] = 1

            if idx == 0:
                current_mass = seq_list[0][0].mass + mass_shift
            else:
                current_mass = current_mass + seq_list[idx][0].mass

            frag_dri = []
            for num in range(0, hn_num + 1):
                new_dict = copy.copy(mod_dict)
                if num != 0:
                    new_dict['HexNAc'] = num
                frag = Fragment(type, idx, new_dict, current_mass)
                frag_dri.append(frag)

            yield frag_dri

    def drop_modification(self, pos, mod_type):
        '''
            Drop a modification by name from a specific residue
        '''
        dropped_index = None
        for i, mod in enumerate(self.seq[pos][1]):
            if mod_type == mod.name:
                dropped_index = i
                break
        try:
            drop_mod = self.seq[pos][1].pop(dropped_index)
            self.mass -= drop_mod.mass
        except:
            raise Exception("Modification not found! %s @ %s" % (mod_type, pos))

    def add_modification(self, pos, mod_type):
        if (pos == -1) or (pos >= len(self.seq)):
            warnings.warn("Invalid modification! Negative Pos: %s Exceeded Length: %s" %
                          ((pos == -1), pos >= len(self.seq)))
            return -1

        mod = Modification(rule=mod_type, mod_pos=pos)
        self.seq[pos][1].append(mod)
        # self.add_mass(pos,mod.mass)
        self.mass += mod.mass

    def append_modification(self, mod):
        if mod.position == -1 | mod.position >= len(self.seq):
            warnings.warn("Invalid modification!")
            return

        pos = mod.position
        self.seq[pos][1].append(mod)
        # self.add_mass(pos,mod.mass)
        self.mass += mod.mass

    def get_sequence(self, start=0):
        """
           Generate user readable sequence string.
           DQYELLC(Carbamidomethyl)LDN(HexNAc)TR
        """

        #seq_str = ''.join([''.join([x.symbol, '(', '|'.join(i.name for i in y), ')']) for x, y in self.seq])
        seq_list = []
        for x, y in self.seq[start:]:
            mod_str = ''
            if y != []:
                mod_str = '|'.join(mod.serialize() for mod in y)
                mod_str = ''.join(['(', mod_str, ')'])
            seq_list.append(''.join([x.symbol, mod_str]))
        return ''.join(seq_list)

    def get_sequence_list(self, start=0):
        return self.seq[start:]

    def at(self, position):
        '''
            A getter retained for backwards compatibility. See __getitem__.
        '''
        return self.seq[position]

    def get_mass(self):
        '''
            A getter retained for backwards compatibility. Defines no special behavior to access mass.
        '''
        return self.mass

    def clone(self):
        '''
            Perform a deep copy of the sequence object by reserializing it and
            then parsing the result.
        '''
        return Sequence(self.get_sequence())
