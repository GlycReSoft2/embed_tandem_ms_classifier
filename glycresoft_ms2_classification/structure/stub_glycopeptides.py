from .modification import Modification
from .composition import Composition
from .sequence import Sequence
from .glycans import oxonium_ions
import re

Deamidation = 0.9840099999999978
Carbamidomethyl = 57.0214
Proton = 1.007276035
HexNAc = 203.07937
Hex = 162.05282
dHex = 146.05791
Water = 18.0105647000


class StubGlycopeptide:

    """Calculates peptide and stub glycopeptide ions, also returns oxonium ions based on glycan composition"""

    def __init__(self, sequence, modification, sites_num, glycan_comp):

        sequence_obj = Sequence(sequence)
        for i in range(0, len(sequence_obj)):
            try:
                sequence_obj.drop_modification(i, "HexNAc")
            except ValueError:
                pass
        sequence = str(sequence_obj)
        self.pep_seq = sequence
        self.mass = sequence_obj.mass
        self.glyco_num = int(sites_num)
        #self.mod = modification
        self.glycan_compo = glycan_comp
        if glycan_comp == '':
            self.dHex = ''
            self.Hex = ''
            self.HexNAc = ''
            self.NeuAc = ''
        else:
            self.glycans = re.findall(r'\b\d+\b', self.glycan_compo)
            self.dHex = int(self.glycans[0])
            self.Hex = int(self.glycans[1])
            self.HexNAc = int(self.glycans[2])
            self.NeuAc = int(self.glycans[3])

        # if modification == '' or modification is None:
        #     self.n_mod = ''
        #     self.mod_type = ''
        # else:
        #     match_pattern = re.compile(r'(\d+)(\w+)')
        #     mods = match_pattern.findall(self.mod)
        #     self.n_mod = int(mods[0][0])
        #     self.mod_type = mods[0][1]

        while True:
            self.begin, self.end = self.pep_seq.find(
                '('), self.pep_seq.find(')')
            if self.begin != -1 and self.end != -1:
                # add this to capture the carbamidomethyls--> + " " +
                # peptide[begin+1:end]
                self.pep_seq = self.pep_seq[
                    :self.begin] + self.pep_seq[self.end + 1:]
            else:
                break

        self.length = len(self.pep_seq)

        # The following part was to add oxonium ions based on presence of sialylation.
        # now we'll instead check if a high mannose composition has 274 or 292 and reduce score:
        # if self.NeuAc == 0:
        #     self.Oxo_ions = [204.0864,186.0754,163.0601,168.0650,138.0542,366.1394]

        # elif self.NeuAc > 0:
        #     self.Oxo_ions = [204.0864,186.0754,163.0601,168.0650,138.0542,366.1394,274.0920,292.1026]

    @classmethod
    def from_sequence(cls, sequence_obj):
        seq = str(sequence_obj)
        num_sites = sequence_obj.mod_index['HexNAc']
        glycan_comp = sequence_obj.glycan
        return StubGlycopeptide(seq, None, num_sites, glycan_comp)

    def get_stubs(self):
        """returns a list of dicts with peptide and stub glycopeptide ions """
        # for i in self.pep_seq:
        #     res = Residue()
        #     res.by_symbol(i)
        #     if i == 'C':
        #         self.mass += Carbamidomethyl
        #     self.mass += res.mass

        # if self.mod_type == "Deamidated":
        #     self.mass += (Deamidation * self.n_mod)
        # elif self.mod_type == '':
        #     pass
        stubs = []

        stubs.append({"key": "peptide", "mass":  (self.mass + Proton)})

        sites = self.glyco_num

        if sites == 0:
            pass
        elif sites > 0:
            if self.dHex == 0:
                for site_i in range(1, sites + 1, 1):
                    for num_hexnac in range(1, 3, 1):
                        if num_hexnac == 1:
                            key = "pep+{num_hexnac}HexNAc+0Hexose-{site_i}sites".format(**locals())
                            # key = "pep+" + \
                            #     str(num_hexnac) + "HexNAc+" + "0Hexose" + "-" + str(
                            #         site_i) + "sites"
                            mass = self.mass + Proton + (site_i * (num_hexnac * HexNAc))
                            stubs.append({"key": key, "mass": mass})

                        elif (num_hexnac == 2):
                            for num_hexose in range(0, 4, 1):
                                key = "pep+{num_hexnac}HexNAc+{num_hexose}Hexose-{site_i}sites".format(**locals())
                                #key = "pep+" + str(num_hexnac) + "HexNAc+" + str(num_hexose) + "Hexose" + "-" + str(site_i) + "sites"
                                mass = self.mass + Proton + (site_i * ((num_hexnac * HexNAc) + (num_hexose * Hex)))
                                stubs.append({"key": key, "mass": mass})
            elif self.dHex > 0:
                for site_i in range(1, sites + 1, 1):
                    for num_hexnac in range(1, 3, 1):
                        key = "pep+{num_hexnac}HexNAc+0Hexose-{site_i}sites".format(**locals())
                        # key = "pep+" + \
                        #     str(num_hexnac) + "HexNAc+" + "0Hexose" + "-" + str(
                        #         site_i) + "sites"
                        mass = self.mass + Proton + (site_i * (num_hexnac * HexNAc))
                        stubs.append({"key": key, "mass": mass})

                        for num_dhex in range(0, 2, 1):
                            for num_hexose in range(0, 4, 1):
                                key = "pep+{num_hexnac}HexNAc+{num_hexose}Hexose\
+{num_dhex}dHex-{site_i}sites".format(**locals())
                                # key = "pep+" + str(num_hexnac) + "HexNAc+" + str(num_hexose) +\
                                #       "Hexose" + str(num_dhex) +\
                                #       "dHex" + "-" + str(site_i) + "sites"
                                mass = self.mass + Proton + \
                                    (site_i * ((num_hexnac * HexNAc) +
                                               (num_hexose * Hex) + (num_dhex * dHex)))
                                stubs.append({"key": key, "mass": mass})

                        # else:
                        #     pass
                        # num_hexnac += 1

        return stubs

    def __iter__(self):
        for stub in self.get_stubs():
            yield stub

    def get_oxonium_ions(self):
        '''returns a list of oxonium ions based on glycan composition'''
        Ox_ions = [{"key": k, "mass": v} for k, v in oxonium_ions.items()]

        return Ox_ions
