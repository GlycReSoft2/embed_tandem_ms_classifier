import copy
import itertools
import random
import re
import time

from math import fabs

from structure.sequence import sequence_to_mass
from structure.sequence import sequence_tokenizer
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


def count_missed_cleavages(sequence, cleavage_end):
    tokens = sequence_tokenizer(sequence)
    missed_cleavages = 0
    for resid, mod in tokens:
        if resid in cleavage_end:
            missed_cleavages += 1
    return missed_cleavages


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


def generate_random_glycopeptides(target_mass, ppm_error=10e-6,
                                  count=20, constant_modifications=None,
                                  variable_modifications=None, glycans=None,
                                  min_length=0, cleavage_start=None, cleavage_end=None,
                                  max_missed_cleavages=0):
    '''
    Given a target mass value and a tolerance threshold around it, create a set of random glycopeptides
    that satisfy the mass requirements.
    '''
    if glycans is None:
        glycans = mammalian_glycans
    if constant_modifications is None:
        constant_modifications = []
    else:
        constant_modifications = copy.deepcopy(constant_modifications)
    if variable_modifications is None:
        variable_modifications = copy.deepcopy(modification.ModificationTable.bootstrap()).values()
    else:
        variable_modifications = copy.deepcopy(variable_modifications)

    if cleavage_start is None or len(cleavage_start) == 0:
        cleavage_start = [""]

    if cleavage_end is None or len(cleavage_end) == 0:
        cleavage_end = [""]

    loc_sequence_to_mass = sequence_to_mass
    loc_fabs = fabs
    loc_gen_glycoforms = generate_random_glycoforms
    joiner = ''.join

    variable_modifications = [mod for mod in variable_modifications if mod.name != "HexNAc"]
    constant_modifications = [mod for mod in constant_modifications if mod.name != "HexNAc"]

    components = generate_component_set(constant_modifications, variable_modifications)
    amino_acid_components = components[:]
    amino_acid_components = set(map(''.join, itertools.product(amino_acid_components + [""], repeat=4)))
    components.extend(generate_n_linked_sequons())
    smallest_mass = min(map(sequence_to_mass, components))

    building_sequences = map(''.join, itertools.product(cleavage_start, components, components))
    random.shuffle(building_sequences)
    accepted_sequences = set()
    accept = accepted_sequences.add
    # Tracks each pass's minimum mass observed. If the minimum mass observed exceeds
    # the maximum acceptable mass, we've run out of solutions.
    min_mass = 0.0
    last_time = start_time = time.time()
    resets = 0
    while resets < 3 and len(accepted_sequences) < count:
        next_time = time.time()
        print(len(accepted_sequences), min_mass, next_time - last_time)
        for seq in accepted_sequences:
            print(seq)
        last_time = next_time

        # If we've not filled the quota, reset and start building a new solution set.
        if loc_fabs((min_mass - target_mass) / float(target_mass)) > ppm_error and min_mass > target_mass:
            building_sequences = map(joiner, itertools.product(cleavage_start, components, components))
            resets += 1
            print("Reset!")

        # Keep the number of random sequences in memory down
        current_min_mass = float('inf')
        try:
            #random.shuffle(building_sequences)
            building_sequences = random.sample(building_sequences, 10000)
        except:
            pass
        random.shuffle(components)
        next_sequences = building_sequences
        print("Combinatorics: %s * %s = %s" %
              (len(building_sequences), len(components), len(building_sequences) * len(components)))

        building_sequences = []
        append_build = building_sequences.append
        extend_build = building_sequences.extend
        # Consider every currently growing sequence combined with every possible amino acid building block
        for sequence in itertools.product(next_sequences, components):
            sequence = joiner(sequence)
            missed_cleavages = count_missed_cleavages(sequence, cleavage_end)
            if missed_cleavages > max_missed_cleavages:
                continue
            sequence_realizations = [sequence] + loc_gen_glycoforms(sequence, glycans)
            for sequence_inst in sequence_realizations:
                try:
                    mass = loc_sequence_to_mass(sequence_inst)
                    #print("Crunching %s (%s)" % (sequence_inst, mass))
                except Exception, e:
                    print("An error occurred during massing of %s" % sequence_inst)
                    print(e)
                    raise
                for sequence_inst_capped in itertools.product([sequence_inst], cleavage_end):
                    inst_mass = mass + loc_sequence_to_mass(sequence_inst_capped[1])
                    # Decide whether the current sequence is a possible solution
                    error = loc_fabs((inst_mass - target_mass) / float(target_mass))
                    if error < ppm_error and is_glycosylated(sequence_inst):
                        accepted = joiner(sequence_inst_capped)
                        print(accepted)
                        accept(accepted)
            # Decide whether the current sequence can be used to build a
            # longer possible solution. Don't save glycoforms.
            if (mass + smallest_mass) < target_mass:
                if(mass < (target_mass/2)):
                    extend_build(joiner(p) for p in itertools.product([sequence_inst], amino_acid_components))
                else:
                    append_build(sequence_inst)

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


# Mass is not being passed quite right yet.
def random_glycopeptide_to_sequence_space(sequence, glycan_string=None, obs_mass=None, ppm_error=0, ms1_score=0):
    seq_obj = Sequence(sequence)
    glycan_map = {}
    modifications = []
    for i, (aa, mods) in enumerate(seq_obj):
        for mod in mods:
            if mod.name in {"Glycan", "HexNAc"}:
                glycan_map[i] = mod.name
            else:
                # Construct the set of acceptable reasons why this modification is here.
                # Infer the least permissive modification rule.
                why = mod.why_valid(aa, i)
                modifications.append(modification.Modification(why, (i,)))

    # Remove glycans from the sequence string to conform to the SequenceSpace expectations
    for site, glycan in glycan_map.items():
        # Don't discard anonymous HexNAcs. Downstream functions can handle them
        if glycan != "HexNAc":
            seq_obj.drop_modification(site, glycan)
    glycosylation_sites = glycan_map.keys()
    if glycan_string is None:
        # Build the semicolon separated string for glycan compositions
        glycan_composition = []
        glycan_composition = [map(int, glycan.replace("Glycan", '').replace("[", "").replace("]", "").split(";"))
                              for glycan in glycan_map.values()]
        glycan_composition = map(sum, zip(*glycan_composition))
        glycan_composition_restring = "[" + ";".join(map(str, glycan_composition)) + "]"
    else:
        glycan_composition_restring = glycan_string
    # Begin generating fragment ions
    b_type = seq_obj.get_fragments('B')
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

    y_type = seq_obj.get_fragments('Y')
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

    peptide_seq = strip_modifications(seq_obj.get_sequence())
    pep_stubs = StubGlycopeptide(peptide_seq, None,
                                 len(glycosylation_sites), glycan_composition_restring)
    stub_ions = pep_stubs.get_stubs()
    oxonium_ions = pep_stubs.get_oxonium_ions()
    ions = {
        "MS1_Score": ms1_score,
        "Obs_Mass": obs_mass if obs_mass is not None else seq_obj.get_mass(),
        "Calc_mass": None,
        "ppm_error": ppm_error,
        "Peptide": peptide_seq,
        "Peptide_mod": None,
        "Glycan": glycan_composition_restring,
        "vol": None,
        "glyco_sites": len(glycan_map),
        "startAA": None,
        "endAA": None,
        "Seq_with_mod": seq_obj.get_sequence(),
        "bare_b_ions": b_ions,
        "b_ions_with_HexNAc": b_ions_HexNAc,
        "bare_y_ions": y_ions,
        "y_ions_with_HexNAc": y_ions_HexNAc,
        "pep_stub_ions": stub_ions,
        "Oxonium_ions": oxonium_ions,
        "Glycopeptide_identifier": seq_obj.get_sequence() + glycan_composition_restring
    }
    return ions


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
    if args.outfile is None:
        args.outfile = args.input_file[:-3] + "random-glycopeptides.fa"
    with open(args.outfile, "w") as fh:
        for mass, count, seqs in sequences:
            for enum, seq in enumerate(seqs):
                fh.write(">" + "%s|%s %s\n" % (mass, count, enum))
                fh.write(seq + '\n')
    print(args.outfile)
