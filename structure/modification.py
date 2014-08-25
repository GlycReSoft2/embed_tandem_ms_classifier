import csv
from copy import deepcopy
import os
import re

from collections import defaultdict

from residue import residue_symbol
# Not used?
from composition import Composition

mod_table = {
    'Carbamidomethyl': ['Cys', 57.0214],
    'Deamidated': ['Asn', 0.9840099999999978],
    'HexNAc': ['Asn', 203.07937],
    'pyroGlu': ['Gln', -17.02655]}


target_string_pattern = re.compile(r'((?P<n_term>N-term)|(?P<c_term>C-term)|(?P<amino_acid>[A-Z]+))')
protein_prospector_name_cleaner = re.compile(r'(?P<name>.*)\s(?P<target>\(.+\))?$')


def extract_targets_from_string(target_string):
    params = target_string_pattern.search(target_string).groupdict()
    amino_acid_targets = params.pop("amino_acid", None)
    amino_acid_targets = map(lambda x: residue_symbol[x], list(amino_acid_targets)) if amino_acid_targets is not None else None
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
        try: # Attempt to parse the protein prospector name which contains the target in
            self.name = protein_prospector_name_cleaner.search(protein_prospector_name).groupdict()['name']
        except:
            self.name = modification_name
        self.unimod_name = modification_name
        self.mass = float(monoisotopic_mass)
        self.protein_prospector_name = protein_prospector_name

    def valid_site(self, amino_acid=None, position_modifiers=None):
        return any([target.valid_site(amino_acid, position_modifiers) for
                    target in self.targets])

    def find_valid_sites(self, sequence):
        valid_indices = []
        position_modifier_rules = {
            0: "N-term",
            (len(sequence) - 1): "C-term"
        }
        for index in range(len(sequence)):
            amino_acid = sequence.at(index)[0].name
            position_modifiers = position_modifier_rules.get(index)
            if(self.valid_site(amino_acid, position_modifiers)):
                valid_indices.append(index)

        return valid_indices

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


class ModificationTable(dict):

    '''Represents the set of modification rules to apply when generating
    the theoretical sequence space.'''

    # Class Constants
    # Path to the bootstrapping data file for Protein Prospector-defined modifications
    _table_definition_file = os.path.join(os.path.dirname(__file__), "data", "ProteinProspectorModifications-for_gly2.csv")

    # Other objects that are described as Modifications that don't appear in Protein Prospector output
    other_modifications = {
        "HexNAc": ModificationRule("N", "HexNAc", "HexNAc", 203.07937),
        "pyroGlu": ModificationRule("Q", "pyroGlu", "pyroGlu", -17.02655)
    }

    # Class Methods
    @classmethod
    def load_from_file(cls, path):
        '''Load the rules definitions from a CSV file and instantiate a ModificationTable from it'''
        definitions = load_modifications_from_file(path)
        return ModificationTable(definitions)

    bootstrapped = None

    @classmethod
    def bootstrap(cls, reuse=True):
        '''Load the bootstrapping file's definitions and instantiate a ModificationTable from it'''
        if cls.bootstrapped is not None and reuse:
            return cls.bootstrapped
        instance = cls.load_from_file(cls._table_definition_file)
        cls.bootstrapped = instance
        return instance

    # Instance Methods
    def __init__(self, modification_definitions=None, constant_modifications=None, variable_modifications=None):
        '''Creates a table from a list of ModificationRules. Assumes that every possibility is available.
        Is an instance of a dictionary object, and therefore can be used as an iterable.

        modification_definitions - A list of ModificationRule objects which define all possible rules from
                                   external sources. They will be filtered down.
        '''
        super(ModificationTable, self).__init__()
        self.modification_rules = []
        for definition in modification_definitions:
            self.modification_rules.append(ModificationRule(**definition))
        self.modification_rules.extend(ModificationTable.other_modifications.values())
        self.constant_modifications = dict()
        self.variable_modifications = dict()
        # A simple look-up by Protein Prospector name to make populating from
        # the GlycReSoft output easier
        self.by_protein_prospector_name = {
            modification.protein_prospector_name: modification for modification in self.modification_rules
        }

        for name, mod in ModificationTable.other_modifications.items():
            self.add_modification(name, mod)
        for mod in self.modification_rules:
            self.add_modification(mod.name, mod)

    def add_modification(self, key, value):
        try:
            mod = self.by_protein_prospector_name[key]
            if mod.name in self:
                self[mod.name] += mod
            else:
                self[mod.name] = mod
        except KeyError:
            try:
                matches = [mod for mod in self.modification_rules if mod.name == key]
                mod = sum(matches[1:], matches[0])
                if mod.name in self:
                    self[mod.name] += mod
                else:
                    self[mod.name] = mod
            except:
                raise

    def get_modification(self, name, position=-1, number=1):
        try:
            mod_rule = self[name]
        except:
            for key, val in self.items():
                print(key, val)
            raise
        modification_instance = Modification(rule=mod_rule, mod_pos=position, mod_num=number)
        return modification_instance


class RestrictedModificationTable(ModificationTable):

    '''Creates a table from a list of ModificationRules. Assumes that we begin with a set of all rules,
    and then selects a subset of them for including in determining valid modifications.

    Is an instance of a dictionary object, and therefore can be used as an iterable.

    modification_definitions - A list of ModificationRule objects which define all possible rules from
                               external sources. They will be filtered down.
    constant_modifications   - A list of modification names present in the table to include in the constant_modifications
                               sub-table
    variable_modifications   - A list of modification names present in the table to include in the variable_modifications
                               sub-table
    '''
    # Class Methods
    @classmethod
    def load_from_file(cls, path, constant_modifications=None, variable_modifications=None):
        '''Load the rules definitions from a CSV file and instantiate a RestrictedModificationTable from it'''
        definitions = load_modifications_from_file(path)
        return RestrictedModificationTable(definitions, constant_modifications, variable_modifications)

    bootstrapped = None

    @classmethod
    def bootstrap(cls, constant_modifications=None, variable_modifications=None, reuse=True):
        '''Load the bootstrapping file's definitions and instantiate a RestrictedModificationTable from it'''
        if cls.bootstrapped is not None and reuse:
            return cls.bootstrapped
        instance = cls.load_from_file(cls._table_definition_file, constant_modifications, variable_modifications)
        cls.bootstrapped = instance
        return instance

    def __init__(self, modification_rules, constant_modifications=None, variable_modifications=None):
        super(RestrictedModificationTable, self).__init__(modification_rules)
        if constant_modifications is not None:
            self._add_constant_modifications(*constant_modifications)
        if variable_modifications is not None:
            self._add_variable_modifications(*variable_modifications)

        # Add the other common modifications as variable modifications so they are in a table
        self._add_variable_modifications(*[x.name for x in ModificationTable.other_modifications.values()])
        self._collect_available_modifications()


    def _add_constant_modifications(self, *protein_prospector_names):
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
                    matches = [mod for mod in self.modification_rules if mod.name == prot_name]
                    mod = sum(matches[1:], matches[0])
                    if mod.name in self.constant_modifications:
                        self.constant_modifications[mod.name] += mod
                    else:
                        self.constant_modifications[mod.name] = mod
                except:
                    pass

    def _add_variable_modifications(self, *protein_prospector_names):
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
                    matches = [mod for mod in self.modification_rules if mod.name == prot_name]
                    mod = sum(matches[1:], matches[0])
                    if mod.name in self.variable_modifications:
                        self.variable_modifications[mod.name] += mod
                    else:
                        self.variable_modifications[mod.name] = mod
                except:
                    pass

    def _collect_available_modifications(self):
        '''Clears the rules stored in the core dictionary (on self, but not its data members),
        and copies all information in the data member sub-tables into the core dictionary'''
        self.clear()
        for name, mod in self.constant_modifications.items():
            self[name] = mod
        for name, mod in self.variable_modifications.items():
            if name in self:
                self[name] += mod
            else:
                self[name] = mod


class Modification:
    """Represents a molecular modification, which may be bound at a given position,
    or may be present at multiple locations. This class apparently pulls double duty,
    and when describing the total number of modifications of a type on a molecule, its
    position attributes are unused."""

    _table = ModificationTable.bootstrap(False)

    def __init__(self, rule, mod_pos=-1, mod_num=1):
        if(isinstance(rule, str)):
            try:
                rule = Modification._table[rule]
            except KeyError:
                raise Exception("Could not resolve Modification Name: %s" % rule)
        self.name = rule.name
        self.mass = rule.mass
        self.position = mod_pos
        self.number = mod_num
        self.rule = rule
        self.target = rule.targets[0].amino_acid_targets[0]

    def valid_site(self, amino_acid=None, position_modifiers=None):
        return self.rule.valid_site(amino_acid, position_modifiers)

    def find_valid_sites(self, sequence):
        return self.rule.find_valid_sites(sequence)

    def __repr__(self):
        rep = "{name}[{number}]@{position}:{rule}".format(**self.__dict__)
        return rep
