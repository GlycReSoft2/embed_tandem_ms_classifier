import copy
import itertools
import random
import re
import time

import multiprocessing

from structure.sequence import sequence_to_mass
from structure.sequence import Sequence
from structure.sequence import strip_modifications
from structure.sequence import sequence_tokenizer

from structure.sequence_space import SequenceSpace
from structure.stub_glycopeptides import StubGlycopeptide

from structure import modification
from structure import residue
from structure import glycan_masses as glycan_masses_provider

glycan_masses = glycan_masses_provider.load_from_file()
min_glycan_mass = min([float(i) for i, j in glycan_masses])

residue_symbols = residue.symbol_to_residue


def generate_n_linked_sequons():
    wild_interval = set(residue_symbols)
    wild_interval.discard("P")
    sequons = map(''.join, itertools.product(["N"], wild_interval, ["S", "T"]))
    return sequons


def generate_component_set(allowed_modifications, modifications_fully_random=False):
    components = list(residue_symbols)
    for mod in allowed_modifications:
        for symbol, name in residue_symbols.items():
            if mod.valid_site(name):
                # If considering all possible modifications, the number of modifications gets
                # a bit ridiculous, so don't include absolutely everything.
                if (modifications_fully_random and random.random() > .2) or (not modifications_fully_random):
                    components.append(symbol + "(" + mod.serialize() + ")")
    return components


def find_n_linked_glycosites(sequence):
    matches = re.finditer(r"[^\)]N[^\(P][ST]", sequence)
    # Glycan is attached to the N, which is the second position in the pattern
    sites = [m.start() + 1 for m in matches]
    return sites


def is_glycosylated(sequence):
    return bool(re.search("Glycan", sequence))


def generate_random_glycoforms(sequence):
    glycoforms = []
    glyco_sites = find_n_linked_glycosites(sequence)
    for n_sites_iter in range(len(glyco_sites)):
        for instance_sites in itertools.combinations(glyco_sites, n_sites_iter):
            glycan = random.choice(glycan_masses)
            glycans_to_add = [modification.AnonymousModificationRule("Glycan" + glycan[1], glycan[0])
                              for s in instance_sites]
            mod_seq = sequence
            for within_site_iter in range(len(instance_sites)):
                pos = instance_sites[within_site_iter] + 1
                mod_seq = mod_seq[:pos] + "(" + glycans_to_add[within_site_iter].serialize() + ")" + mod_seq[pos:]
            glycoforms.append(mod_seq)
    return glycoforms


def evaluate_sequence(sequence, target_mass, lower_bound, upper_bound):
    current_min_mass = float('inf')
    accepted_sequences = []
    building_sequences = []
    sequence = ''.join(sequence)
    sequence_realizations = [sequence] + generate_random_glycoforms(sequence)
    for sequence_inst in sequence_realizations:
        mass = sequence_to_mass(sequence_inst)
        sequence_string = sequence_inst

        # Decide whether the current sequence is a possible solution
        if ((target_mass - lower_bound) < mass < (target_mass + upper_bound)) and\
           is_glycosylated(sequence_string):
            accepted_sequences.append(sequence_string)

        # Decide whether the current sequence can be used to build a
        # longer possible solution
        if mass < (target_mass + upper_bound):
            building_sequences.append(sequence_string)

        # If this sequence has a lower mass than the smallest mass observed,
        # update the smallest mass observed for this pass
        if mass < current_min_mass:
            current_min_mass = mass

    return accepted_sequences, building_sequences, current_min_mass


def generate_random_glycopeptides(target_mass, lower_bound, upper_bound,
                                  count=20, allowed_modifications=None, glycan_map=None):
    '''
    Given a target mass value and a tolerance threshold around it, create a set of random glycopeptides
    that satisfy the mass requirements.
    '''
    if glycan_map is None:
        glycan_map = dict()
    modifications_fully_random = False
    if allowed_modifications is None:
        modifications_fully_random = True
        allowed_modifications = copy.deepcopy(modification.ModificationTable.bootstrap()).values()
    else:
        allowed_modifications = copy.deepcopy(allowed_modifications)
    allowed_modifications = [mod for mod in allowed_modifications if mod.name != "HexNAc"]
    components = generate_component_set(allowed_modifications, modifications_fully_random)
    components.extend(generate_n_linked_sequons())
    #building_sequences = components[:]
    building_sequences = map(''.join, itertools.product(components[:], repeat=2))

    accepted_sequences = []

    # Tracks each pass's minimum mass observed. If the minimum mass observed exceeds
    # the maximum acceptable mass, we've run out of solutions.
    min_mass = 0.0
    last_time = start_time = time.time()
    resets = 0
    while resets < 3 and len(accepted_sequences) < count:
        next_time = time.time()
        print(len(accepted_sequences), min_mass, next_time - last_time)
        last_time = next_time

        # If we've not filled the quota, reset and start building a new solution set.
        if min_mass > (target_mass + upper_bound):
            building_sequences = map(''.join, itertools.product(components[:], repeat=2))
            resets += 1
        # Keep the number of random sequences in memory down
        current_min_mass = float('inf')
        try:
            random.shuffle(building_sequences)
            building_sequences = random.sample(building_sequences, 5000)
        except:
            pass
        random.shuffle(components)
        next_sequences = building_sequences
        building_sequences = []
        for sequence in itertools.product(next_sequences, components):
            sequence = ''.join(sequence)
            sequence_realizations = [sequence] + generate_random_glycoforms(sequence)
            for sequence_inst in sequence_realizations:
                mass = sequence_to_mass(sequence_inst)
                sequence_string = sequence_inst

                # Decide whether the current sequence is a possible solution
                if ((target_mass - lower_bound) < mass < (target_mass + upper_bound)) and\
                   is_glycosylated(sequence_string):
                    accepted_sequences.append(sequence_string)

                # Decide whether the current sequence can be used to build a
                # longer possible solution
                if mass < (target_mass + upper_bound):
                    building_sequences.append(sequence_string)

                # If this sequence has a lower mass than the smallest mass observed,
                # update the smallest mass observed for this pass
                if mass < current_min_mass:
                    current_min_mass = mass

            # Update the tracking of the smallest mass observed
            min_mass = current_min_mass

    print(len(accepted_sequences))
    if len(accepted_sequences) > count:
        accepted_sequences = random.sample(accepted_sequences, count)
    print(time.time() - start_time)
    return accepted_sequences


def random_glycopeptide_to_sequence_space(sequence):
    seq_obj = Sequence(sequence)
    glycan_map = {}
    modifications = []
    for i, (aa, mods) in enumerate(seq_obj):
        for mod in mods:
            if "Glycan" in mod.name:
                glycan_map[i] = mod.name
            else:
                # Construct the set of acceptable reasons why this modification is here.
                # Infer the least permissive modification rule.
                why = mod.why_valid(aa, i)
                modifications.append(modification.Modification(why, (i,)))

    # Remove glycans from the sequence string to conform to the SequenceSpace expectations
    for site, glycan in glycan_map.items():
        seq_obj.drop_modification(site, glycan)
    glycosylation_sites = glycan_map.keys()
    # Build the semicolon separated string for glycan compositions
    glycan_composition = []
    glycan_composition = [map(int, glycan.replace("Glycan", '').replace("[", "").replace("]", "").split(";"))
                          for glycan in glycan_map.values()]
    glycan_composition = map(sum, zip(*glycan_composition))
    glycan_composition_restring = "[" + ";".join(map(str, glycan_composition)) + "]"

    sequence_space = SequenceSpace(seq_obj.get_sequence(), None, glycosylation_sites, modifications)
    results = []
    for seq_inst in sequence_space.get_theoretical_sequence(len(glycosylation_sites)):
        b_type = seq_inst.get_fragments('B')
        b_ions = []
        b_ions_HexNAc = []
        for b in b_type:
            for fm in b:
                key = fm.get_fragment_name()
                if key == "B1" or re.search(r'B1\+', key):
                    # B1 Ions aren't actually seen in reality, but are an artefact of the generation process
                    # so do not include them in the output
                    continue
                mass = fm.get_mass()
                if "HexNAc" in key:
                    b_ions_HexNAc.append({"key": key, "mass": mass})
                else:
                    b_ions.append({"key": key, "mass": mass})

        y_type = seq_inst.get_fragments('Y')
        y_ions = []
        y_ions_HexNAc = []
        for y in y_type:
            for fm in y:
                key = fm.get_fragment_name()
                mass = fm.get_mass()
                if "HexNAc" in key:
                    y_ions_HexNAc.append({"key": key, "mass": mass})
                else:
                    y_ions.append({"key": key, "mass": mass})

        pep_stubs = StubGlycopeptide(strip_modifications(seq_inst.get_sequence()), None,
                                     len(glycosylation_sites), glycan_composition_restring)
        stub_ions = pep_stubs.get_stubs()
        oxonium_ions = pep_stubs.get_oxonium_ions()
        ions = {
            "precursor_mass": seq_inst.get_mass(),
            "bare_b_ions": b_ions,
            "b_ions_with_HexNAc": b_ions_HexNAc,
            "bare_y_ions": y_ions,
            "y_ions_with_HexNAc": y_ions_HexNAc,
            "pep_stub_ions": stub_ions,
            # Why this is capitalized, nobody knows.
            "Oxonium_ions": oxonium_ions
        }
        results.append(ions)
    return results
