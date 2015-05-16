import itertools
import random

from bitarray import bitarray

from glycresoft_ms2_classification.ms import spectra
from glycresoft_ms2_classification.structure import sequence, glycans, stub_glycopeptides


oxonium_ions = dict(glycans.oxonium_ions)
extract_annotations = spectra.extract_annotations
Sequence = sequence.Sequence
StubGlycopeptide = stub_glycopeptides.StubGlycopeptide


def sequence_to_bit_vector(seq_obj, frag_type='b', tol=1):
    tol = int(tol)
    total_mass = int(seq_obj.mass * tol)
    bit_vector = bitarray('0') * (total_mass)
    for fragment in itertools.chain.from_iterable(seq_obj.get_fragments(frag_type)):
        bit_vector[int(fragment.mass) * tol] = True
    return bit_vector


def spectrum_to_bit_vector(spec_obj, tol=1):
    tol = int(tol)
    total_mz = max(
        int(max(spec_obj.tandem_data, key=lambda x: x.neutral_mass).neutral_mass) + 1,
        int(spec_obj.neutral_mass) + 1)
    bit_vector = bitarray('0') * total_mz * tol
    for tand in spec_obj.tandem_data:
        bit_vector[int(tand.neutral_mass) * tol] = True
    mask_oxonium_ions(bit_vector, tol)
    return bit_vector


def sequences_from_annotation(precursor):
    seqs = map(Sequence, {s for a in extract_annotations(*precursor.tandem_data)
                          for s, k in a.items()})
    return ((precursor, s) for s in seqs)


def mask_oxonium_ions(precursor_bits, tol=1):
    for ion_type, mass in oxonium_ions.items():
        precursor_bits[int(mass * tol)] = 0


def mask_masses(mass_list, precursor_bits, tol=1):
    for mass in mass_list:
        precursor_bits[int(mass * tol)] = 0


def preprocess(psms, frag_type='b', tol=1):
    for spec, pep in psms:
        stubs = [f['mass'] for f in StubGlycopeptide.from_sequence(pep).get_stubs()]
        pep_bits = sequence_to_bit_vector(pep, frag_type, tol)
        spec_bits = spectrum_to_bit_vector(spec, tol)
        mask_masses(stubs, spec_bits, tol)
        yield spec_bits, pep_bits


def run(psms, frag_type='b', tol=1):
    psms_bits = list(preprocess(
        itertools.chain.from_iterable(map(sequences_from_annotation, psms)),
        frag_type, tol
        ))
    alpha = find_alpha(psms_bits)
    beta = find_beta(psms_bits)
    prior_fragmentation = find_empirical_random(psms_bits)
    return alpha, beta, prior_fragmentation


def truncate(spec, pep):
    '''Discard bits from spec exceeding pep's length'''
    size = len(pep)
    return spec[:size]


def find_alpha(psms):
    matched = 0
    total = 0
    for spec, pep in psms:
        spec = truncate(spec, pep)
        paired = (spec & pep).count()
        if paired == 0:
            continue
        matched += paired
        total += pep.count()
    return matched / float(total)


def find_beta(psms):
    unmatched = 0
    potential = 0
    for spec, pep in psms:
        spec = truncate(spec, pep)
        paired = (spec & pep).count()
        if paired == 0:
            continue
        unmatched += spec.count() - paired
        # Assumes that the total mass of the sequence will be spanned by its fragments
        # This is untrue for glycopeptide precursors.
        potential += len(pep) - pep.count()
    return unmatched / float(potential)


def find_empirical_random(psms):
    count = 0
    total = 0
    for spec, pep in psms:
        count += random.choice(pep)
        total += 1
    return count / float(total)
