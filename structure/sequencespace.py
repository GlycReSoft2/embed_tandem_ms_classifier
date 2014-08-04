from sequence import Sequence
from operator import and_
from functools import reduce
from modification import Modification
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
        self.seq = Sequence(seq) # Sequence object
        self.glycan_composition = glycan_compo
        self.candidate_sites = glycan_sites
        self.modifications = mod_list

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

        ## Get the candidate sites for all modification
        for mod in self.modifications:
            if mod.position != -1: # The position specified.
                ix_bound.append([mod.position]) # Single element list
            elif mod.target!= '': # The target specified.
                ix_list = [ix for ix in range(self.seq.length) if self.seq.at(ix)[0].name == mod.target]
                ## temp_list has format like [(1,2,3), (2,3,4)]
                temp_list = [ix for ix in itertools.combinations(ix_list, mod.number)]
                ix_bound.append(temp_list)
            else:
                raise Exception('Unqualified modification!')

        ## Initialize the choice index for each modification type.
        indices = [0] * n

        while True:
            if n != 0:
                for i in reversed(range(n)):
                    ## If not achiving the last choice of current index
                    if indices[i] != len(ix_bound[i]): # Within boundary, just out of the loop
                        break
                    else: # Out of boundary, reset the index.
                        indices[i] = 0
                        if i > 0:
                            indices[i-1] += 1
                else:
                    return seq_space

            ## Check if current indecies are qualifed.
                ix_sites = [ix_bound[ss][indices[ss]] for ss in range(n)]
            else:
                ix_sites = []

            common_sites = set().union(*ix_sites)
            glyco_sites = set(self.candidate_sites).difference(common_sites)
            #glyco_num = glyco_compo['HexNAc']
            if len(common_sites) != sum(map(len,ix_sites)) or (num_sites > len(glyco_sites)): # Invalid config.
                indices[i] += 1
                continue

            raw_seq = copy.deepcopy(self.seq)

            for x in range(n):
                for mod_site in ix_bound[x][indices[x]]:
                    raw_seq.addModification(mod_site, self.modifications[x].name)

            ## Get available glycosylation sites.
            
            #upper_limit = (min(max_sites, len(glyco_sites)) if max_sites > 0 else len(glyco_sites))

            #for m in range(1, upper_limit+1):
            for sites in itertools.combinations(glyco_sites, num_sites):
                temp_seq = copy.deepcopy(raw_seq)
                # Append HexNAc to the corresponding sites.
                for site in sites:
                    gly_mod = Modification("HexNAc", site, 1, Residue("HexNAc").mass, 'Asn')
                    temp_seq.appendModification(gly_mod)
                seq_space.append(temp_seq)

            if n == 0:
                return seq_space

            # Only increase the last index.
            indices[-1] += 1




