import copy
import itertools

from structure.modification import Modification, ModificationTable
from structure.sequence import Sequence


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
modifications:{modifications}""".format(seq2=self.seq.get_sequence(), **self.__dict__)
        return rep

    def get_modification_sites(self):
        modification_index_bound = []
        for mod in self.modifications:
            if mod.position != -1:
                # The modification position is pre-specified
                modification_index_bound.append([mod.position])
            else:
                valid_sites = mod.find_valid_sites(self.seq)
                all_combinations = [position for position in itertools.combinations(valid_sites, mod.number)]
                if(len(all_combinations) == 0):
                    raise NoSitesFoundException("No Valid Sites Found: %r %r in\n %r" % (mod, valid_sites, self.seq))
                modification_index_bound.append(all_combinations)

        return modification_index_bound

    def compose_sequence(self, modification_sites, num_sites):

        seq_space = []
        num_modifications = len(self.modifications)
        modification_indices = [0] * num_modifications
        while(True):
            if num_modifications != 0:
                for ind in reversed(range(num_modifications)):
                    if(modification_indices[ind] != len(modification_sites[ind])):
                        break
                    else:
                        modification_indices[ind] = 0
                        if(ind > 0):
                            modification_indices[ind - 1] += 1
                # There are no modification index combinations, just return the sequence.
                else:
                    return seq_space

                combination_index_sites = [modification_sites[site][modification_indices[site]] for
                                           site in range(num_modifications)]
            else:
                combination_index_sites = []
            common_sites = set(combination_index_sites)
            glycosylation_sites = set(self.candidate_sites).difference(common_sites)

            if(len(common_sites) != sum(map(len, combination_index_sites)) or
               (num_sites > len(glycosylation_sites))):  # Invalid Configuration
                modification_indices[ind] += 1
                continue
                raise Exception("Invalid Configuration", modification_sites,
                                modification_indices,
                                combination_index_sites,
                                common_sites)

            raw_sequence = copy.deepcopy(self.seq)

            for index in range(num_modifications):
                for mod_site in modification_sites[index][modification_indices[index]]:
                    raw_sequence.add_modification(mod_site, self.modifications[index].rule)

            for sites in itertools.combinations(glycosylation_sites, num_sites):
                temp_seq = copy.deepcopy(raw_sequence)
                for site in sites:
                    gly_mod = Modification(ModificationTable.other_modifications["HexNAc"], site, 1)
                    temp_seq.append_modification(gly_mod)
                seq_space.append(temp_seq)

            if num_modifications == 0:
                return seq_space

            modification_indices[-1] += 1

    def get_theoretical_sequence(self, num_sites):
        modification_index_bound = self.get_modification_sites()
        sequence_space = self.compose_sequence(modification_index_bound, num_sites)
        return sequence_space


class NoSitesFoundException(Exception):
    pass


class UnqualifiedModifierException(Exception):
    pass
