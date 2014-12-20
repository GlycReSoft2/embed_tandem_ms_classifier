from .modification import Modification


class Fragment:
    """ Fragment from glycopeptide """
    def __init__(self, type, pos, mod_dict, mass):
        self.type = type
        # The mass value is the bare backbone's mass
        self.mass = mass

        self.pos = pos
        self.mod_dict = mod_dict

        for key, value in self.mod_dict.items():
            self.mass += Modification(key).mass * value

    def get(self):
        """
            simply return string like B2, Y3 with no modificaiton information.
        """
        fragment_name = []
        fragment_name.append(self.type)
        fragment_name.append(str(self.pos))
        return ''.join(fragment_name)

    def get_mass(self):
        return self.mass

    def get_modification_number(self, mod_name):
        if mod_name in self.mod_dict:
            return self.mod_dict[mod_name]
        else:
            return 0

    def get_modifications(self):
        for name, number in self.mod_dict.items():
            yield [name, number]

    def get_fragment_name(self):
        """
           Connect the information into string.
        """
        fragment_name = []
        fragment_name.append(self.type)
        fragment_name.append(str(self.pos))

        ## Only concerned modification is reported.
        concerned_mod = ['HexNAc']
        for mod_name in concerned_mod:
            if mod_name in self.mod_dict:
                if self.mod_dict[mod_name] > 1:
                    fragment_name.extend(['+', str(self.mod_dict[mod_name]), mod_name])
                elif self.mod_dict[mod_name] == 1:
                    fragment_name.extend(['+', mod_name])
                else:
                    pass

        return ''.join(fragment_name)

    def __repr__(self):
        return "Fragment(%(type)s @ %(pos)s %(mass)s [%(mod_dict)s])" % self.__dict__
