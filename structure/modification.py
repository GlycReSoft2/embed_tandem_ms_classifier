import csv
from copy import deepcopy
import re

from collections import defaultdict

# Not used?
from composition import Composition

mod_table = {
    'Carbamidomethyl': ['Cys', 57.0214],
    'Deamidated': ['Asn', 0.9840099999999978],
    'HexNAc': ['Asn', 203.07937],
    'pyroGlu': ['Gln', -17.02655]}


target_string_pattern = re.compile(r'((?P<n_term>N-term)|(?P<c_term>C-term)|(?P<amino_acid>[A-Z]+))')


def extract_targets_from_string(target_string):
    params = target_string_pattern.search(target_string).groupdict()
    amino_acid_targets = params.pop("amino_acid", None)
    amino_acid_targets = list(amino_acid_targets) if amino_acid_targets is not None else None
    position_modifiers = [pos for pos, flag in params.items() if flag is not None] if len(params) != 0 else None
    position_modifiers = position_modifiers if len(position_modifiers) != 0 else None

    target = ModificationTarget(amino_acid_targets, position_modifiers)
    return target


class ModificationTarget(object):

    '''Specifies the range of targets that a given modification can be found at.'''

    def __init__(self, site_targets=None, position_modifiers=None):
        self.amino_acid_targets = site_targets
        self.position_modifiers = position_modifiers

    def valid_site(self, amino_acid=None, position_modifiers=None):
        valid = False
        # Validate amino acid target target
        valid = (self.amino_acid_targets is None) or (amino_acid in self.amino_acid_targets)
        valid = valid and ((self.position_modifiers is None) or
                          (position_modifiers in self.position_modifiers))

        return valid

    def __repr__(self):
        rep = "{amino_acid_targets}@{position_modifiers}".format(**self.__dict__)
        return rep


class ModificationRule(object):

    '''Represent the name, mass, and position specifities associated with a given modification.
    Additionally stores information on the represented modification in Protein Prospector.'''

    def __init__(self, amino_acid_specificity, modification_name,
                 protein_prospector_name, monoisotopic_mass, **kwargs):
        self.targets = [extract_targets_from_string(amino_acid_specificity)]
        self.name = modification_name
        self.mass = monoisotopic_mass
        self.protein_prospector_name = protein_prospector_name

    def valid_site(self, amino_acid=None, position_modifiers=None):
        return any([target.valid_site(amino_acid, position_modifiers) for
                    target in self.targets])

    def find_valid_sites(self, peptide_fragment):
        pass

    def __repr__(self):
        rep = "{name}({protein_prospector_name}):{mass}:@{targets}".format(**self.__dict__)
        return rep

    def __add__(self, other):
        if(isinstance(other, ModificationRule)):
            dup = deepcopy(self)
            dup.targets.extend(other.targets)
            return dup
        else:
            raise TypeError("Unsupported types. Can only add two ModificationRules together.")


def load_modifications_from_file(path):
    modification_definitions = []
    with open(path, 'rb') as fh:
        parser = csv.DictReader(fh)
        modification_definitions = [rec for rec in parser]
    return modification_definitions


class ModificationTable(object):

    @classmethod
    def load_from_file(cls, path):
        definitions = load_modifications_from_file(path)
        return ModificationTable(definitions)

    def __init__(self, modification_definitions, constant_modifications=None, variable_modifications=None):
        self.all_modifications = []
        for definition in modification_definitions:
            if("Protein prospector names" in definition):
                definition['protein_prospector_name'] = definition.pop("Protein prospector names")
            self.all_modifications.append(ModificationRule(**definition))
        self.constant_modifications = dict()
        self.variable_modifications = dict()

        self.by_protein_prospector_name = {
            modification.protein_prospector_name: modification for modification in self.all_modifications
        }
        if constant_modifications is not None:
            self.add_constant_modifications(constant_modifications)
        if variable_modifications is not None:
            self.add_variable_modifications(variable_modifications)

    def add_constant_modifications(self, *protein_prospector_names):
        for prot_name in protein_prospector_names:
            try:
                mod = self.by_protein_prospector_name[prot_name]
                if mod.name in self.constant_modifications:
                    self.constant_modifications[mod.name] += mod
                else:
                    self.constant_modifications[mod.name] = mod
            except KeyError, e:
                print(e)
                try:
                    matches = [mod for mod in self.all_modifications if mod.name == prot_name]
                    mod = sum(matches[1:], matches[0])
                    if mod.name in self.constant_modifications:
                        self.constant_modifications[mod.name] += mod
                    else:
                        self.constant_modifications[mod.name] = mod
                except:
                    pass

    def add_variable_modifications(self, *protein_prospector_names):
        for prot_name in protein_prospector_names:
            try:
                mod = self.by_protein_prospector_name[prot_name]
                if mod.name in self.variable_modifications:
                    self.variable_modifications[mod.name].targets.extend(mod.targets)
                else:
                    self.variable_modifications[mod.name] = mod
            except KeyError, e:
                print(e)
                try:
                    matches = [mod for mod in self.all_modifications if mod.name == prot_name]
                    mod = sum(matches[1:], matches[0])
                    if mod.name in self.variable_modifications:
                        self.variable_modifications[mod.name] += mod
                    else:
                        self.variable_modifications[mod.name] = mod
                except:
                    pass


class Modification:

    """description of class"""

    def __init__(self, mod_name, mod_pos=-1, mod_num=1, mass=0.0, target=''):
        self.name = mod_name
        self.position = mod_pos
        self.number = mod_num
        key = mod_table.get(mod_name)
        if key is None:
            self.mass = mass
            self.target = target
            mod_table[mod_name] = [target, mass]
        else:
            self.mass = mod_table.get(mod_name)[1]
            self.target = mod_table.get(mod_name)[0]
