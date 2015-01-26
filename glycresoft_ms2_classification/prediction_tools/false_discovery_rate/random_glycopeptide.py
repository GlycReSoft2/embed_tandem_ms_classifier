import copy
import itertools
import random
import re
import logging
import multiprocessing
logger = logging.getLogger("random_sequence_builder")
from math import fabs


from glycresoft_ms2_classification.structure.sequence import Sequence
from glycresoft_ms2_classification.structure.sequence import GrowingSequence
from glycresoft_ms2_classification.structure.sequence import Protease
from glycresoft_ms2_classification.structure.sequence import strip_modifications
from glycresoft_ms2_classification.structure.sequence import sequence_to_mass
from glycresoft_ms2_classification.structure.stub_glycopeptides import StubGlycopeptide
from glycresoft_ms2_classification.structure.composition import Composition
from glycresoft_ms2_classification.structure import modification
from glycresoft_ms2_classification.structure import residue
from glycresoft_ms2_classification.structure import glycans as glycan_lib
from glycresoft_ms2_classification.utils.mass_heap import MassHeap

mammalian_glycans = glycan_lib.load_from_file()

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
            if mod.valid_site(name, "internal"):
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


def has_glycan(sequence):
    for mod_name, count in sequence.mod_index.items():
        if ("Glycan" in mod_name) or ("HexNAc" in mod_name) and (count > 0):
            return count
    return False


def generate_random_glycopeptides(target_mass, ppm_error=10e-6, count=20, constant_modifications=None,
                                  variable_modifications=None, glycans=None, min_length=0, cleavage_start=None,
                                  cleavage_end=None, max_missed_cleavages=1, max_glycosylations=2):
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
        variable_modifications = []
    else:
        variable_modifications = copy.deepcopy(variable_modifications)

    if cleavage_start is None or len(cleavage_start) == 0:
        cleavage_start = [""]

    if cleavage_end is None or len(cleavage_end) == 0:
        cleavage_end = [""]

    cleavage_pattern = Protease(cleavage_start, cleavage_end)

    variable_modifications = [
        mod for mod in variable_modifications if mod.name != "HexNAc"]
    constant_modifications = [
        mod for mod in constant_modifications if mod.name != "HexNAc"]

    components = MassHeap(map(lambda x: GrowingSequence(x, cleavage_pattern), generate_component_set(
        constant_modifications, variable_modifications)))
    sequons = MassHeap(
        map(lambda x: GrowingSequence(x, cleavage_pattern),
            itertools.chain.from_iterable(
                map(lambda x: ("{0}({1}){2}".format(x[0], g.as_modification().serialize(), x[1:])
                               for g in glycans),
                    generate_n_linked_sequons()
                    )
            )
            )
    )

    loc_fabs = fabs
    water = Composition("H2O").mass

    def reset_target_mass():
        return (water + target_mass) - min(p.mass for p in candidate.pad())
    solutions = set()
    max_iter = count * 10000
    iter_count = 0
    candidate = GrowingSequence("", cleavage_pattern)
    mass_to_meet = reset_target_mass()
    while(len(solutions) < count and iter_count < max_iter):
        can_glycosylate = (len(candidate) > min_length / 3) and \
            (has_glycan(candidate) < max_glycosylations) and \
            (random.random() > .7)
        options = list(components.get_lower_than(mass_to_meet))

        if(can_glycosylate):
            glycosylated_options = list(sequons.get_lower_than(mass_to_meet))
            options += glycosylated_options

        #logger.debug("%s options for extension, mass to meet: %s, %s" % (len(options), mass_to_meet, str(candidate)))
        next_part = random.choice(options)
        candidate.extend(next_part)
        mass_to_meet -= (next_part.mass - water)
        # print(str(candidate), candidate.missed_cleavages, len(candidate))
        # Reset, too many missed cleavages?
        if candidate.missed_cleavages > max_missed_cleavages:
            #print("Too many missed cleavages: %s, Reset!" % candidate.missed_cleavages)
            candidate = GrowingSequence("", cleavage_pattern)
            mass_to_meet = reset_target_mass()

        for padded_sequence in candidate.pad():
            # Only consider glycosylated sequences
            if has_glycan(candidate) < 1:
                break

            # Only consider longer sequences
            if(len(padded_sequence) < min_length):
                continue

            error = loc_fabs(
                (target_mass - padded_sequence.mass) / float(target_mass))
            # logger.debug("%s, %s, %s" %
            #              (padded_sequence, padded_sequence.mass, error))
            # Accept?
            if error <= ppm_error:
                #logger.debug("Accepting %s %s" %
                #             (padded_sequence, padded_sequence.mass))
                solutions.add(str(padded_sequence))

        # Reset, too big?
        if mass_to_meet < components[0].mass:
            candidate = GrowingSequence("", cleavage_pattern)
            mass_to_meet = reset_target_mass()

        iter_count += 1

    return solutions


def forge_prediction_record(sequence, reference):
    seq_obj = Sequence(sequence) if isinstance(
        sequence, basestring) else copy.deepcopy(sequence)

    calculated_mass = seq_obj.mass

    glycan_map = {}
    modifications = []
    loc_rules = modification.get_position_modifier_rules_dict(seq_obj)
    for i, (aa, mods) in enumerate(seq_obj):
        for mod in mods:
            if "Glycan" in mod.name:
                glycan_map[i] = mod.name
            else:
                # Construct the set of acceptable reasons why this modification is here.
                # Infer the least permissive modification rule.
                try:
                    why = mod.why_valid(aa, loc_rules[i])
                    modifications.append(modification.Modification(why, (i,)))
                except AttributeError:
                    print(mod)
                    raise

    # Remove glycans from the sequence string to conform to the SequenceSpace
    # expectations
    for site, glycan in glycan_map.items():
        # Don't discard anonymous HexNAcs. Downstream functions can handle them
        if glycan != "HexNAc":
            seq_obj.drop_modification(site, glycan)
            seq_obj.add_modification(site, "HexNAc")

    # Build the semicolon separated string for glycan compositions
    glycan_composition = []
    glycan_composition = [map(int, glycan.replace("Glycan", '').replace("[", "").replace("]", "").split(";"))
                          for glycan in glycan_map.values()]
    glycan_composition = map(sum, zip(*glycan_composition))
    glycan_composition_restring = "[" + \
        ";".join(map(str, glycan_composition)) + "]"

    forgery = reference.copy()
    forgery["_old_Glycopeptide_identifier"] = str(sequence)
    forgery.Calc_mass = calculated_mass
    forgery.Obs_Mass = forgery.Calc_mass - reference.ppm_error
    forgery.Glycopeptide_identifier = str(seq_obj)
    forgery.Glycan = glycan_composition_restring
    forgery.glyco_sites = len(glycan_map)
    forgery.Seq_with_mod = str(seq_obj)

    return forgery


# Mass is not being passed quite right yet.
def random_glycopeptide_to_sequence_space(sequence, proxy, glycan_string=None):
    seq_obj = Sequence(sequence) if isinstance(
        sequence, basestring) else sequence
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

    # Remove glycans from the sequence string to conform to the SequenceSpace
    # expectations
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
        glycan_composition_restring = "[" + \
            ";".join(map(str, glycan_composition)) + "]"
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
        "MS1_Score": proxy["MS1_Score"],
        "Obs_Mass": seq_obj.mass - proxy["ppm_error"],
        "Calc_mass": seq_obj.mass,
        "ppm_error": proxy["ppm_error"],
        "Peptide": peptide_seq,
        "Peptide_mod": None,
        "Glycan": glycan_composition_restring,
        "vol": proxy['vol'],
        "glyco_sites": len(glycan_map),
        "startAA": -1,
        "endAA": -1,
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


class RandomGlycopeptideBuilder(object):

    def __init__(self, ppm_error, constant_modifications=None, variable_modifications=None,
                 glycans=None, cleavage_start=None, cleavage_end=None, ignore_sequences=None):
        if glycans is None:
            glycans = mammalian_glycans
        self.glycans = set(glycans)
        if constant_modifications is None:
            constant_modifications = []
        else:
            constant_modifications = copy.deepcopy(constant_modifications)
        self.constant_modifications = constant_modifications
        if variable_modifications is None:
            variable_modifications = []
        else:
            variable_modifications = copy.deepcopy(variable_modifications)
        self.variable_modifications = variable_modifications

        if cleavage_start is None or len(cleavage_start) == 0:
            cleavage_start = [""]
        self.cleavage_start = cleavage_start

        if cleavage_end is None or len(cleavage_end) == 0:
            cleavage_end = [""]
        if ignore_sequences is None:
            ignore_sequences = []
        self.ignore_sequences = set(ignore_sequences)
        self.cleavage_end = cleavage_end

        self.ppm_error = ppm_error

        self.cleavage_pattern = Protease(cleavage_start, cleavage_end)

        self.variable_modifications = [
            mod for mod in variable_modifications if mod.name != "HexNAc"]
        self.constant_modifications = [
            mod for mod in constant_modifications if mod.name != "HexNAc"]

        self.components = MassHeap(map(lambda x: GrowingSequence(x, self.cleavage_pattern), generate_component_set(
            constant_modifications, variable_modifications)))
        self.sequons = MassHeap(
            map(lambda x: GrowingSequence(x, self.cleavage_pattern),
                itertools.chain.from_iterable(
                    map(lambda x: ("{0}({1}){2}".format(x[0], g.as_modification().serialize(), x[1:])
                                   for g in self.glycans),
                        generate_n_linked_sequons()))
                )
        )

    def generate_random(self, target_mass, count, max_missed_cleavages=3, max_glycosylations=2, min_length=5):
        loc_fabs = fabs
        water = Composition("H2O").mass

        components = self.components
        sequons = self.sequons
        cleavage_pattern = self.cleavage_pattern
        ppm_error = self.ppm_error

        def reset_target_mass():
            return (water + target_mass) - min(p.mass for p in candidate.pad())
        solutions = set()
        max_iter = count * 10000
        iter_count = 0
        candidate = GrowingSequence("", self.cleavage_pattern)
        mass_to_meet = reset_target_mass()
        while(len(solutions) < count and iter_count < max_iter):
            try:
                can_glycosylate = (len(candidate) > min_length / 3) and \
                    (has_glycan(candidate) < max_glycosylations) and \
                    (random.random() > .7)
                options = list(components.get_lower_than(mass_to_meet))

                if(can_glycosylate):
                    glycosylated_options = list(
                        sequons.get_lower_than(mass_to_meet))
                    options += glycosylated_options
                if len(options) == 0:
                    candidate = GrowingSequence("", cleavage_pattern)
                    mass_to_meet = reset_target_mass()

                next_part = random.choice(options)
                candidate.extend(next_part)
                mass_to_meet -= (next_part.mass - water)
                # print(str(candidate), candidate.missed_cleavages, len(candidate))
                # Reset, too many missed cleavages?
                if candidate.missed_cleavages > max_missed_cleavages:
                    #print("Too many missed cleavages: %s, Reset!" % candidate.missed_cleavages)
                    candidate = GrowingSequence("", cleavage_pattern)
                    mass_to_meet = reset_target_mass()

                for padded_sequence in candidate.pad():
                    # Only consider glycosylated sequences
                    if has_glycan(candidate) < 1:
                        break

                    # Only consider longer sequences
                    if(len(padded_sequence) < min_length):
                        continue

                    error = loc_fabs(
                        (target_mass - padded_sequence.mass) / float(target_mass))
                    #logger.debug("%s, %s, %s" %
                    #             (padded_sequence, padded_sequence.mass, error))
                    # Accept?
                    if error <= ppm_error:
                        if str(padded_sequence) in self.ignore_sequences:
                            continue
                        #logger.debug("Accepting %s %s" %
                                    #(padded_sequence, padded_sequence.mass))
                        solutions.add(str(padded_sequence))

                # Reset, too big?
                if mass_to_meet < components[0].mass:
                    candidate = GrowingSequence("", cleavage_pattern)
                    mass_to_meet = reset_target_mass()

                iter_count += 1
            except IndexError, e:
                pname = multiprocessing.current_process().name
                logger.error("[%s] RandomGlycopeptideBuilder: Exception occurred", pname, exc_info=e)
        return solutions
