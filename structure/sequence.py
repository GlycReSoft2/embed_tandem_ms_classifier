from modification import Modification
from composition import Composition
from residue import Residue
from fragment import Fragment

import copy
import re
import warnings


class Sequence:

    """Name and mass of residues (glycan/peptide)"""

    def __init__(self, sequence):
        # a. parse the sequence.
        match_pattern = re.compile(r'([A-Z])(?:\((\w+)\))?')  # Use non-capturing group
        seq_list = match_pattern.findall(sequence)
        self.mass = 0.0
        self.seq = []
        for item in seq_list:
            res = Residue()
            res.bySymbol(item[0])
            if item[1] != '':
                mod = Modification(item[1])
        # Update!!! The mod mass should not be integrated into residue mass.
        #res.mass += mod.mass
                self.seq.append([res, [mod]])
        # Update!!! The mod mass is appended to the precursor directly.
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

    def addMass(self, pos, mass):
        """
                This function will direcly modify the mass of the residue instead of appending the modificaion.
                Use it with caution!
            """
        self.seq[pos][0].mass += mass
        self.mass += mass

    def getFragments(self, type):
        """ Return a list of mass values"""

        mass_shift = 0.0

        # HexNAc number
        #hn_num = 0

        # The set of modification names.
        mod_dict = {}

        # The key is the position, the value is an array of fragments.
        # And the first element is always bare fragment.
        # The total number of HexNAc on the fragment should be recorded.\

        #fragment_dict = {}
        fragment_series = []
        if type == 'B':
            #mass_shift = Composition('H').mass - Composition('e').mass
            seq_list = self.seq

            mass_shift = Composition('H').mass - Composition('e').mass
            #fragment_series.append(seq_list[0][0].mass + Composition('H').mass - Composition('e').mass)

        elif type == 'Y':
            mass_shift = Composition('H2O').mass + Composition('H').mass - Composition('e').mass
            #mass_shift = Composition('OH').mass - Composition('e').mass
            seq_list = list(reversed(self.seq))
            #fragment_series.append(seq_list[0][0].mass + Composition('H2O').mass + Composition('H').mass - Composition('e').mass)

        #bare_frag_mass = seq_list[0][0].mass + mass_shift

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
            #frag = Fragment(type, idx, mod_dict, current_mass)
            # frag_dri.append(frag)
            for num in range(0, hn_num + 1):
                # It's OK if the HexNAc is 0.
                new_dict = copy.copy(mod_dict)
                if num != 0:
                    new_dict['HexNAc'] = num
                frag = Fragment(type, idx, new_dict, current_mass)
                frag_dri.append(frag)

            # fragment_series.append(frag_dri)
            yield frag_dri

        # for item in seq_list[1:-1]:
        #    res_mass = item[0].mass
        #    fragment_series.append(fragment_series[-1] + res_mass)

        # return fragment_series

    def addModification(self, pos, mod_type):
        if (pos == -1) or (pos >= len(self.seq)):
            warnings.warn("Invalid modification!")
            return

        mod = Modification(rule=mod_type, mod_pos=pos)
        self.seq[pos][1].append(mod)
        # self.addMass(pos,mod.mass)
        self.mass += mod.mass

    def appendModification(self, mod):
        if mod.position == -1 | mod.position >= len(self.seq):
            warnings.warn("Invalid modification!")
            return

        pos = mod.position
        self.seq[pos][1].append(mod)
        # self.addMass(pos,mod.mass)
        self.mass += mod.mass

    def getSequence(self):
        """
           Generate user readable sequence string.
           DQYELLC(Carbamidomethyl)LDN(HexNAc)TR
        """

        #seq_str = ''.join([''.join([x.symbol, '(', '|'.join(i.name for i in y), ')']) for x, y in self.seq])
        seq_list = []
        for x, y in self.seq:
            mod_str = ''
            if y != []:
                mod_str = '|'.join(mod.name for mod in y)
                mod_str = ''.join(['(', mod_str, ')'])
            seq_list.append(''.join([x.symbol, mod_str]))
        return ''.join(seq_list)

    def getSequenceList(self):
        return self.seq

    def at(self, position):
        return self.seq[position]

    def getMass(self):
        return self.mass
