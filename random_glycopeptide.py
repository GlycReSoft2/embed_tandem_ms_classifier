import copy
import itertools
import random
import re
import time

from structure.sequence import sequence_to_mass
from structure.sequence import Sequence
from structure.sequence import strip_modifications

from structure.sequence_space import SequenceSpace
from structure.stub_glycopeptides import StubGlycopeptide

from structure import modification
from structure import residue
from structure import glycan_masses as glycan_masses_provider

mammalian_glycans = glycan_masses_provider.load_from_file()

residue_symbols = residue.symbol_to_residue
smallest_mass = sequence_to_mass("G")


def generate_n_linked_sequons():
    wild_interval = set(residue_symbols)
    wild_interval.discard("P")
    sequons = map(''.join, itertools.product(["N"], wild_interval, ["S", "T"]))
    return sequons


def generate_component_set(constant_modifications, variable_modifications):
    components = set(residue_symbols)
    const_modified = set()
    for mod in constant_modifications:
        for symbol, name in residue_symbols.items():
            if mod.valid_site(name):
                components.add(symbol + "(" + mod.serialize() + ")")
                const_modified.add(symbol)

    # After applying constant modifications, we can't have the same residue unmodified, so
    # remove it from the component set.
    map(components.discard, const_modified)

    for mod in variable_modifications:
        for symbol, name in residue_symbols.items():
            if mod.valid_site(name):
                components.add(symbol + "(" + mod.serialize() + ")")

    components = list(components)
    return components


def find_n_linked_glycosites(sequence):
    matches = re.finditer(r"N[^\(P][ST]", sequence)
    # Glycan is attached to the N, which is the second position in the pattern
    sites = [m.start() for m in matches]
    return sites


def is_glycosylated(sequence):
    return bool(re.search("Glycan", sequence))


def partition(sequence, indices, offset=0):
    if len(indices) == 0:
        return [sequence]
    segments = []
    current_segment = ""
    seq_iter = 0
    segment_iter = 0
    current_split = indices[segment_iter]
    while(seq_iter < len(sequence) and segment_iter < len(indices)):
        char = sequence[seq_iter]
        current_segment += char
        if seq_iter >= current_split + offset:
            segments.append(current_segment)
            current_segment = ""
            segment_iter += 1
            if segment_iter < len(indices):
                current_split = indices[segment_iter]
                assert(current_split != len(sequence))
            else:
                segments.append(sequence[seq_iter + 1:])
        seq_iter += 1
    return segments


def generate_random_glycoforms(sequence, glycans):
    glycoforms = []
    glyco_sites = find_n_linked_glycosites(sequence)
    for n_sites_iter in range(len(glyco_sites)):
        for instance_sites in itertools.combinations(glyco_sites, n_sites_iter):
            glycan = random.choice(glycans)
            glycans_to_add = [glycan.as_modification()
                              for s in instance_sites]
            mod_seq = ""
            segments = partition(sequence, instance_sites)
            for within_site_iter in range(len(instance_sites)):
                glycan_str = "(%s)" % glycans_to_add[within_site_iter].serialize()
                mod_seq += segments[within_site_iter] + glycan_str
            mod_seq += segments[-1]
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
                                  count=20, constant_modifications=None,
                                  variable_modifications=None, glycans=None):
    '''
    Given a target mass value and a tolerance threshold around it, create a set of random glycopeptides
    that satisfy the mass requirements.
    '''
    if glycans is None:
        glycans = mammalian_glycans
    if constant_modifications is None:
        constant_modifications = copy.deepcopy(modification.ModificationTable.bootstrap()).values()
    else:
        constant_modifications = copy.deepcopy(constant_modifications)

    if variable_modifications is None:
        variable_modifications = copy.deepcopy(modification.ModificationTable.bootstrap()).values()
    else:
        variable_modifications = copy.deepcopy(variable_modifications)

    variable_modifications = [mod for mod in variable_modifications if mod.name != "HexNAc"]
    constant_modifications = [mod for mod in constant_modifications if mod.name != "HexNAc"]

    components = generate_component_set(constant_modifications, variable_modifications)
    components.extend(generate_n_linked_sequons())
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
            sequence_realizations = [sequence] + generate_random_glycoforms(sequence, glycans)
            for sequence_inst in sequence_realizations:
                try:
                    mass = sequence_to_mass(sequence_inst)
                except Exception, e:
                    print("An error occurred during massing of %s" % sequence_inst)
                    print(e)
                    raise
                sequence_string = sequence_inst

                # Decide whether the current sequence is a possible solution
                if ((target_mass - lower_bound) < mass < (target_mass + upper_bound)) and\
                   is_glycosylated(sequence_string):
                    accepted_sequences.append(sequence_string)

                # Decide whether the current sequence can be used to build a
                # longer possible solution
                if mass < (target_mass + upper_bound - smallest_mass):
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


def generate_random_glycopeptides2(target_mass, lower_bound, upper_bound, count, components):
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
            sequence_inst = ''.join(sequence)
            try:
                mass = sequence_to_mass(sequence_inst)
            except Exception, e:
                print("An error occurred during massing of %s" % sequence_inst)
                print(e)
                raise
            # Decide whether the current sequence is a possible solution
            if ((target_mass - lower_bound) < mass < (target_mass + upper_bound)) and\
               is_glycosylated(sequence_inst):
                accepted_sequences.append(sequence_inst)

            # Decide whether the current sequence can be used to build a
            # longer possible solution
            if mass < (target_mass + upper_bound - smallest_mass):
                building_sequences.append(sequence_inst)

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


def build_single_random_glycopeptide(target_mass, lower_bound, upper_bound, constant_modifications=None,
                                     variable_modifications=None, glycans=None):
    sequence = ""
    mass = 0
    if glycans is None:
        glycans = random.sample(mammalian_glycans, 25)
    if constant_modifications is None:
        constant_modifications = copy.deepcopy(modification.ModificationTable.bootstrap()).values()
    else:
        constant_modifications = copy.deepcopy(constant_modifications)

    if variable_modifications is None:
        variable_modifications = copy.deepcopy(modification.ModificationTable.bootstrap()).values()
    else:
        variable_modifications = copy.deepcopy(variable_modifications)

    variable_modifications = [mod for mod in variable_modifications if mod.name != "HexNAc"]
    constant_modifications = [mod for mod in constant_modifications if mod.name != "HexNAc"]

    components = generate_component_set(constant_modifications, variable_modifications)
    components.extend("%s(%s)%s" % (site[:1],
                                    glycan.as_modification().serialize(),
                                    site[1:]) for
                      site, glycan in itertools.product(generate_n_linked_sequons(), glycans))
    component_masses = map(sequence_to_mass, components)
    last_add_mass = 0
    last_piece_len = 0
    sequence = random.choice(components)
    mass = sequence_to_mass(sequence)
    loop_count = 0
    while (mass <= (target_mass + upper_bound)) and loop_count < 100000:
        if(target_mass + upper_bound - mass < smallest_mass):
            mass -= last_add_mass
            sequence = sequence[:-last_piece_len]
        loop_count += 1
        remaining_mass = (target_mass + upper_bound) - mass
        valid_components = [components[i] for i in xrange(len(components)) if component_masses[i] <= remaining_mass]
        next_piece = random.choice(valid_components)
        next_mass = sequence_to_mass(next_piece)
        if (target_mass - lower_bound) <= mass <= (target_mass + upper_bound):
            mass += next_mass
            sequence += next_piece
            break
        mass += next_mass
        last_add_mass = next_mass
        sequence += next_piece
        last_piece_len = len(next_piece)
    if (target_mass - lower_bound) <= mass <= (target_mass + upper_bound) and is_glycosylated(sequence):
        return sequence
    else:
        print("No valid sequence generated")
        return None


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

        peptide_seq = strip_modifications(seq_inst.get_sequence())
        pep_stubs = StubGlycopeptide(peptide_seq, None,
                                     len(glycosylation_sites), glycan_composition_restring)
        stub_ions = pep_stubs.get_stubs()
        oxonium_ions = pep_stubs.get_oxonium_ions()
        ions = {
            "MS1_Score": None,
            "Obs_Mass": seq_inst.get_mass(),
            "Calc_mass": None,
            "ppm_error": None,
            "Peptide": peptide_seq,
            "Peptide_mod": None,
            "Glycan": glycan_composition_restring,
            "vol": None,
            "glyco_sites": len(glycan_map),
            "startAA": None,
            "endAA": None,
            "Seq_with_mod": seq_inst.get_sequence(),
            "bare_b_ions": b_ions,
            "b_ions_with_HexNAc": b_ions_HexNAc,
            "bare_y_ions": y_ions,
            "y_ions_with_HexNAc": y_ions_HexNAc,
            "pep_stub_ions": stub_ions,
            "Oxonium_ions": oxonium_ions,
            "Glycopeptide_identifier": seq_inst.get_sequence() + glycan_composition_restring
        }
        results.append(ions)
    return results


def process_task(param, tolerance=5.0, number=20, *args, **kwargs):
    try:
        mass, count = param
        sequences = generate_random_glycopeptides(mass, tolerance, tolerance, count=count * number)
    except Exception, e:
        sequences = []
        print(e, mass, count)
    return (mass, count, sequences)


def tester(mass, tolerance=1, *args, **kwargs):
    print(mass, tolerance, args, kwargs)

if __name__ == '__main__':
    import argparse
    import os
    from collections import Counter
    from multiprocessing import cpu_count, Pool
    from functools import partial
    app = argparse.ArgumentParser()
    app.add_argument("input_file")
    app.add_argument("-t", "--tolerance", type=float, default=5.0,
                     help="Tolerance range around the parameter mass to accept random glycopeptides")
    app.add_argument("-n", "--number", type=int, default=20, help="Number of random glycopeptides per match")
    app.add_argument("-o", "--outfile", type=str, default=None, required=False)
    app.add_argument("-c", "--cores", type=int, default=1, required=False, help="Number of cores to use")

    args = app.parse_args()
    args.cores = cpu_count() if args.cores > cpu_count() else args.cores
    mass_list = []
    if os.path.splitext(args.input_file)[1] == ".csv":
        import csv
        with open(args.input_file, 'rb') as fh:
            reader = csv.DictReader(fh)
            mass_list.extend(map(float, [r["Obs_Mass"] for r in reader]))
    else:
        with open(args.input_file, 'rb') as fh:
            mass_list.extend(map(float, [r for r in fh]))

    mass_list = map(lambda x: round(x, 4), mass_list)
    quantities = Counter(mass_list)
    random_sequences = {}
    workers = Pool(args.cores)
    task = partial(process_task, tolerance=args.tolerance, number=args.number)
    sequences = workers.map(task, quantities.iteritems(), 25)
    print(len(sequences))
    import pickle
    pickle.dump(sequences, open("backup.pickle", 'wb'))

    if args.outfile is None:
        args.outfile = args.input_file[:-3] + "random-glycopeptides.fa"
    with open(args.outfile, "w") as fh:
        for mass, count, seqs in sequences:
            for enum, seq in enumerate(seqs):
                fh.write(">" + "%s|%s %s\n" % (mass, count, enum))
                fh.write(seq + '\n')
    print(args.outfile)
