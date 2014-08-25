from sequence import Sequence
from operator import and_
from functools import reduce
from modification import Modification, ModificationTable
from residue import Residue

import copy
import itertools
import warnings


class SequenceSpace:

    """Generate all theoretical glycopeptide sequences"""

    def __init__(self, seq, glycan_compo, glycan_sites, mod_list):
        """
            seq -- sequence code
            glycan_compo -- glycan compositions, dict.
            glycan_sites -- sets of candidate sites for glycosylation
            mod_list -- list of modifications.
        """
        # Filter the glycan composition. Get the max number of HexNAc
        self.seq = Sequence(seq)  # Sequence object
        self.glycan_composition = glycan_compo
        self.candidate_sites = glycan_sites
        self.modifications = mod_list

    def __repr__(self):
        rep = """
sequence:{seq2}
glycan:{glycan_composition}
candidate_sites:{candidate_sites}
modifications:{modifications}""".format(seq2=self.seq.getSequence(), **self.__dict__)
        return rep

    def get_modification_sites(self):

        n = len(self.modifications)
        modification_index_bound = []
        for mod in self.modifications:
            if mod.position != -1:
                # The modification position is pre-specified
                modification_index_bound.append([mod.position])
            else:
                valid_sites = mod.find_valid_sites(self.seq)
                all_combinations = [position for position in itertools.combinations(valid_sites, mod.number)]
                if(len(all_combinations) == 0):
                    raise Exception("No Valid Sites Found: %r %r in\n %r" % (mod, valid_sites, self.seq))
                modification_index_bound.append(all_combinations)

        return modification_index_bound

    def compose_sequence(self, modification_sites, num_sites):

        seq_space = []

        num_modifications = len(self.modifications)
        modification_indices = [0] * num_modifications

        while(True):
            if num_modifications != 0:
                #print("num_modifications is not 0")
                for ind in reversed(range(num_modifications)):
                    #print(modification_indices)
                    #print(modification_indices[ind], len(modification_sites[ind]))
                    if(modification_indices[ind] != len(modification_sites[ind])):
                        #print("Breaking")
                        break
                    else:
                        modification_indices[ind] = 0
                        if(ind > 0):
                            modification_indices[ind - 1] += 1
                else:
                    #print("returning 1")
                    #print(seq_space)
                    return seq_space

                combination_index_sites = [modification_sites[site][modification_indices[site]] for site in range(num_modifications)]

            else:
                combination_index_sites = []


            common_sites = set(*combination_index_sites)
            glycosylation_sites = set(self.candidate_sites).difference(common_sites)


            #print("Indices",modification_indices)
            #print("Bounds",modification_sites)
            #print("Sites",combination_index_sites)

            #print("Configuration: ", len(common_sites), sum(map(len, combination_index_sites)),
            #  num_sites,  len(glycosylation_sites))
            if(len(common_sites) != sum(map(len, combination_index_sites)) or
               (num_sites > len(glycosylation_sites))):  # Invalid Configuration
                modification_indices[ind] += 1
                continue
                raise Exception("Invalid Configuration", modification_sites,
                                modification_indices,
                                combination_index_sites,
                                common_sites)

            #print("generating modifications")
            raw_sequence = copy.deepcopy(self.seq)

            for index in range(num_modifications):
                for mod_site in modification_sites[index][modification_indices[index]]:
                    #print(mod_site, self.modifications[index].rule)
                    raw_sequence.addModification(mod_site, self.modifications[index].rule)

            for sites in itertools.combinations(glycosylation_sites, num_sites):
                temp_seq = copy.deepcopy(raw_sequence)
                for site in sites:
                    gly_mod = Modification(ModificationTable.other_modifications["HexNAc"], site, 1)
                    temp_seq.appendModification(gly_mod)
                seq_space.append(temp_seq)
                #print(seq_space)


            if num_modifications == 0:
                #print("returning 2")
                return seq_space

            modification_indices[-1] += 1


    def get_theoretical_sequence(self, num_sites):
        modification_index_bound = self.get_modification_sites()
        #print(modification_index_bound, self.modifications)
        sequence_space = self.compose_sequence(modification_index_bound, num_sites)
        return sequence_space

    def getTheoreticalSequence(self, num_sites):
        """
            Get theoretical sequence tailored for fragmenation
            max_sites -- the number of maximum glycolsylation sites.
            -1 means unlimited.
        """

        #raw_seq = self.seq

        seq_space = []
        occupied_sites = []

        #exploreSequence(mod_set, 0, raw_seq, occupied_sites, seq_space)

        n = len(self.modifications)

        ix_bound = []

        # Get the candidate sites for all modification
        for mod in self.modifications:
            if mod.position != -1:  # The position specified.
                ix_bound.append([mod.position])  # Single element list
            elif mod.target != '':  # The target specified.
                ix_list = [ix for ix in range(self.seq.length) if self.seq.at(ix)[0].name == mod.target]
                # temp_list has format like [(1,2,3), (2,3,4)]
                temp_list = [ix for ix in itertools.combinations(ix_list, mod.number)]
                ix_bound.append(temp_list)
            else:
                raise Exception('Unqualified modification!')

        # Initialize the choice index for each modification type.
        indices = [0] * n

        while True:
            if n != 0:
                #print("n is not 0")
                for i in reversed(range(n)): # i is dumped into the outer scope and referenced later
                    # If not achiving the last choice of current index
                    #print(indices)
                    #print(indices[i], len(ix_bound[i]))
                    if indices[i] != len(ix_bound[i]):  # Within boundary, just out of the loop
                        #print("Breaking")
                        break
                    else:  # Out of boundary, reset the index.
                        indices[i] = 0
                        if i > 0:
                            indices[i - 1] += 1
                else:
                    #print(seq_space)
                    return seq_space

            # Check if current indecies are qualifed.
                ix_sites = [ix_bound[ss][indices[ss]] for ss in range(n)]
            else:
                ix_sites = []

            #print("Indices",indices)
            #print("Bounds",ix_bound)
            #print("Sites",ix_sites)
            common_sites = set().union(*ix_sites)
            glyco_sites = set(self.candidate_sites).difference(common_sites)
            #glyco_num = glyco_compo['HexNAc']

            #print("Configuration: ", len(common_sites), sum(map(len, ix_sites)), num_sites, len(glyco_sites))
            if len(common_sites) != sum(map(len, ix_sites)) or (num_sites > len(glyco_sites)):  # Invalid config.
                indices[i] += 1
                continue

            #print("generating modifications")
            raw_seq = copy.deepcopy(self.seq)

            for x in range(n):
                for mod_site in ix_bound[x][indices[x]]:
                    #print(mod_site, self.modifications[x].rule)
                    raw_seq.addModification(mod_site, self.modifications[x].rule)

            # Get available glycosylation sites.

            #upper_limit = (min(max_sites, len(glyco_sites)) if max_sites > 0 else len(glyco_sites))

            # for m in range(1, upper_limit+1):
            for sites in itertools.combinations(glyco_sites, num_sites):
                temp_seq = copy.deepcopy(raw_seq)
                # Append HexNAc to the corresponding sites.
                for site in sites:
                    #gly_mod = Modification("HexNAc", site, 1, Residue("HexNAc").mass, 'Asn')
                    gly_mod = Modification(ModificationTable.other_modifications["HexNAc"], site, 1)
                    temp_seq.appendModification(gly_mod)
                seq_space.append(temp_seq)
                #print(seq_space)

            if n == 0:
                #print("returning 2")
                return seq_space

            # Only increase the last index.
            indices[-1] += 1
