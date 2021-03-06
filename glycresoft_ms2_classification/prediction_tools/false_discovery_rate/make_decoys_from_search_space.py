import re
import multiprocessing
import functools
import itertools
import logging
import time
import argparse

import sqlitedict

from glycresoft_ms2_classification.structure import sequence, constants
from glycresoft_ms2_classification.structure.parser import sequence_tokenizer_respect_sequons

Sequence = sequence.Sequence
list_to_sequence = sequence.list_to_sequence
logger = logging.getLogger()


def fragments(seq):
    seq_mod = seq.get_sequence()
    fragments = zip(*map(seq.break_at, range(1, len(seq))))
    b_type = fragments[0]
    b_ions = []
    b_ions_HexNAc = []
    for b in b_type:
        for fm in b:
            key = fm.get_fragment_name()
            if key == ("B1" or re.search(r'B1\+', key)) and constants.EXCLUDE_B1:
                # B1 Ions aren't actually seen in reality, but are an artefact of the generation process
                # so do not include them in the output
                continue
            mass = fm.get_mass()
            golden_pairs = fm.golden_pairs
            if "HexNAc" in key:
                b_ions_HexNAc.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})
            else:
                b_ions.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})

    y_type = fragments[1]  # seq.get_fragments('Y')
    y_ions = []
    y_ions_HexNAc = []
    for y in y_type:
        for fm in y:
            key = fm.get_fragment_name()
            mass = fm.get_mass()
            golden_pairs = fm.golden_pairs
            if "HexNAc" in key:
                y_ions_HexNAc.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})
            else:
                y_ions.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})

    fragments_dict = {
        "Seq_with_mod": seq_mod,
        "Glycopeptide_identifier": str(seq),
        "bare_b_ions": b_ions,
        "bare_y_ions": y_ions,
        "b_ions_with_HexNAc": b_ions_HexNAc,
        "y_ions_with_HexNAc": y_ions_HexNAc
    }
    return fragments_dict


def pair_rotate(sequence):
    chunks = []
    gen = iter(sequence)
    while True:
        chunk = []
        try:
            chunk = [gen.next()]
            chunk.append(gen.next())
            chunks.append(chunk)
        except StopIteration:
            if len(chunk) > 0:
                chunks.append(chunk)
            break
    rev_seq = []
    for chunk in reversed(chunks):
        rev_seq.extend(chunk)
    return rev_seq


def make_decoy(theoretical_sequence, prefix_len=0, suffix_len=0):
    seq = sequence_tokenizer_respect_sequons(theoretical_sequence["Glycopeptide_identifier"])
    pref = seq[:prefix_len]
    if suffix_len == 0:
        suf = ""
        body = seq[prefix_len:]
    else:
        suf = seq[-suffix_len:]
        body = seq[prefix_len:-suffix_len]
    body = body[::-1]
    rev_seq = (list_to_sequence(pref + list(body) + suf))
    if str(list_to_sequence(seq)) == str(rev_seq):
        rot_body = pair_rotate(body)
        rev_seq = (list_to_sequence(pref + list(rot_body) + suf))

    fragments_dict = fragments(rev_seq)
    decoy = dict(theoretical_sequence)
    decoy.update(fragments_dict)
    return decoy


def taskmain(predictions_path, n_processes=4, prefix_len=0, suffix_len=0, out=None):
    decoy_path = (predictions_path if out is None else out).rsplit("scored", 1)[0] + "decoy.db"

    predictions_metadata = sqlitedict.open(predictions_path, "metadata")

    metadata = dict(predictions_metadata.items())
    metadata["tag"] = (metadata["tag"] if metadata["tag"] not in ["", None] else "") + "decoy"
    metadata["decoy_ratio"] = 1
    metadata_table = sqlitedict.open(decoy_path, 'metadata')
    metadata_table.update(metadata)
    metadata_table.commit()
    metadata_table.close()

    logger.info("Reading sequence space from %r", predictions_path)
    theoretical_sequences = sqlitedict.open(predictions_path, "theoretical_search_space")

    decoy_table = sqlitedict.open(decoy_path, "theoretical_search_space", journal_mode="OFF")
    logger.info("Writing decoys to %r", decoy_path)
    cntr = 0
    task_fn = functools.partial(make_decoy, prefix_len=prefix_len, suffix_len=suffix_len)
    if n_processes > 1:
        worker_pool = multiprocessing.Pool(n_processes)
        for decoy in itertools.chain(
                worker_pool.imap_unordered(task_fn, theoretical_sequences.itervalues(), chunksize=500)):
            decoy_table[cntr] = decoy
            cntr += 1
            if cntr % 10000 == 0:
                logger.info("Processed %d decoys", cntr)
                decoy_table.commit()
        worker_pool.terminate()
        worker_pool.close()
    else:
        for decoy in itertools.imap(task_fn, theoretical_sequences.itervalues()):
            decoy_table[cntr] = decoy
            cntr += 1
            if cntr % 10000 == 0:
                logger.info("Processed %d decoys", cntr)
                decoy_table.commit()
    decoy_table.commit()
    print len(decoy_table), len(theoretical_sequences)
    if len(decoy_table) == 0:
        raise Exception("No decoys generated")
    logger.info("Decoy creation complete.")
    return decoy_path
