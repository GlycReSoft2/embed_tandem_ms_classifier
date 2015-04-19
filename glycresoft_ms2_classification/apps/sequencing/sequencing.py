import itertools
import logging
from math import fabs
from collections import deque
from functools import partial

from disk_queue import SqliteDiskQueue

from glycresoft_ms2_classification.structure.sequence import sequence_to_mass
from glycresoft_ms2_classification.structure import residue, modification
from glycresoft_ms2_classification.structure.composition import Composition
from glycresoft_ms2_classification.utils.mass_heap import MassHeap
from glycresoft_ms2_classification.structure.sequence import SimplePeptide

residue_symbols = dict(residue.symbol_to_residue)
residue_symbols.pop("I")
logging.basicConfig(level='DEBUG', format='%(asctime)s %(message)s', datefmt='%I:%M:%S %p')
logger = logging.getLogger()


precursor_mass_shift = Composition('H2O').mass
y_mass_shift = Composition('H2O').mass + Composition('H').mass - Composition('e').mass
b_mass_shift = Composition('H').mass - Composition('e').mass


def y_mass(seq):
    return sequence_to_mass(seq) - precursor_mass_shift + y_mass_shift


class SimpleFragment(SimplePeptide):
    def __init__(self, sequence_str, mass=None, missed_cleavages=0, cleaver=None, mod_index=None, length=None):
        if mass is None:
            mass = sequence_to_mass(str(sequence_str)) - precursor_mass_shift
        super(SimpleFragment, self).__init__(
            sequence_str, mass, missed_cleavages, cleaver, mod_index, length)


class SequenceRecord(object):
    def __init__(self, sequence, total_gaps=0, current_gaps=0, index=0, kind="X", graph=None, matches=None):
        self.sequence = sequence
        self.total_gaps = total_gaps
        self.current_gaps = current_gaps
        self.mass = self.sequence.mass
        self.index = index
        self.graph = graph
        self.matches = matches or []
        self.kind = kind

    def extend(self, seq, gap=False):
        return SequenceRecord(
            SimplePeptide(
                str(self.sequence) + str(seq),
                self.mass + seq.mass,
                length=self.sequence.length + seq.length),
            self.total_gaps + gap,
            self.current_gaps + 1 if gap else 0,
            graph=self.graph,
            index=self.index + 1)

    def __repr__(self):
        rep = "<{kind}:{sequence} {mass} {current_gaps}/{total_gaps}>".format(**self.__dict__)
        return rep


class ResultsGroup(object):
    def __init__(self, sequences, parts):
        self.sequences = sorted(sorted(sequences, key=lambda x: -x.mass), key=lambda x: x.total_gaps)
        self.parts = parts


def ppm_error(query, match):
    return (query - match)/match


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


def comparator_loop(seq, upper_threshold, lower_threshold, parts, tandem, max_running_gaps, max_total_gaps):
    match = []
    for msms in reversed(list(tandem.get_higher_than(lower_threshold))):
        if msms.neutral_mass > upper_threshold:
            break
        for part in parts:
            mass_query = part.mass + y_mass_shift + seq.mass
            if seq.kind != 'b' and (fabs(ppm_error(msms.neutral_mass, mass_query)) < 2e-5):
                ext = seq.extend(part, False)
                ext.matches.append(msms)
                ext.kind = 'y'
                yield (ext)
                match.append(msms)
                # logger.info("Match on %r -> %r", ext, msms)
            mass_query += -y_mass_shift + b_mass_shift
            if seq.kind != 'y' and (fabs(ppm_error(msms.neutral_mass, mass_query)) < 2e-5):
                ext = seq.extend(part, False)
                ext.kind = 'b'
                yield ext
                match.append(msms)
                # logger.info("Match on %r -> %r", ext, msms)
    if len(match) == 0:
        # logger.info("No match on %r", seq)
        if seq.current_gaps + 1 <= max_running_gaps and seq.total_gaps + 1 <= max_total_gaps:
            for part in parts:
                yield (seq.extend(part, True))


def sequence_spectra(ms_spectrum, drop_mass=0, constant_modifications=None, variable_modifications=None,
                     max_running_gaps=1, max_total_gaps=2):
    constant_modifications = constant_modifications or []
    variable_modifications = variable_modifications or []
    variable_modifications += [modification.Modification("HexNAc").rule]
    precursor_mass = ms_spectrum.neutral_mass - drop_mass
    logger.info("Precursor Mass: %f", precursor_mass)
    previous_sequences = SqliteDiskQueue()
    parts = map(SimpleFragment,
                generate_component_set(constant_modifications,
                                       variable_modifications))
    tandem = MassHeap(ms_spectrum.tandem_data)

    # Get starting component
    match = False
    for part in parts:
        for msms in tandem:
            if (fabs(ppm_error(msms.neutral_mass, part.mass + y_mass_shift)) < 2e-5):
                previous_sequences.append(SequenceRecord(part, kind='y'))
                match = True
            if (fabs(ppm_error(msms.neutral_mass, part.mass + b_mass_shift)) < 2e-5):
                previous_sequences.append(SequenceRecord(part, kind='b'))
                match = True
    if not match:
        for part in parts:
            previous_sequences.append(SequenceRecord(part, 1, 1))
    next_sequences = SqliteDiskQueue()
    solutions = deque(maxlen=4 * max_total_gaps)
    min_part = min(parts, key=lambda x: x.mass).mass
    max_part = max(parts, key=lambda x: x.mass).mass
    while len(previous_sequences) > 0:
        for seq in previous_sequences:
            match = []
            lower = (seq.mass + min_part)
            upper = seq.mass + max_part
            lower_threshold = lower - (y_mass_shift + lower * 2e-5)
            upper_threshold = upper + (y_mass_shift + upper * 2e-5)
            for msms in reversed(list(tandem.get_higher_than(lower_threshold))):
                if msms.neutral_mass > upper_threshold:
                    break
                for part in parts:
                    mass_query = part.mass + y_mass_shift + seq.mass
                    if seq.kind != 'b' and (fabs(ppm_error(msms.neutral_mass, mass_query)) < 2e-5):
                        ext = seq.extend(part, False)
                        ext.matches.append(msms)
                        ext.kind = 'y'
                        next_sequences.append(ext)
                        match.append(msms)
                        # logger.info("Match on %r -> %r", ext, msms)
                    mass_query += -y_mass_shift + b_mass_shift
                    if seq.kind != 'y' and (fabs(ppm_error(msms.neutral_mass, mass_query)) < 2e-5):
                        ext = seq.extend(part, False)
                        ext.kind = 'b'
                        next_sequences.append(ext)
                        match.append(msms)
                        # logger.info("Match on %r -> %r", ext, msms)
            if len(match) == 0:
                # logger.info("No match on %r", seq)
                if seq.current_gaps + 1 <= max_running_gaps and seq.total_gaps + 1 <= max_total_gaps:
                    for part in parts:
                        next_sequences.append(seq.extend(part, True))
        logger.info("Round over, %d candidates", len(next_sequences))
        if len(next_sequences) == 0:
            return ResultsGroup([seq for round in solutions for seq in round], parts)
        previous_sequences = next_sequences
        solutions.append(seq for seq in next_sequences if len(seq.matches) > 0)
        next_sequences = SqliteDiskQueue()
