import logging
decoy_logger = logging.getLogger("make_decoys")
import argparse
import multiprocessing
import functools
import itertools
import math
import sqlitedict
import random
import re
import copy
import time

from glycresoft_ms2_classification.structure import modification, sequence, stub_glycopeptides
from glycresoft_ms2_classification.structure.glycans import from_predictions as glycans_from_predictions
from glycresoft_ms2_classification import prediction_tools
from .random_glycopeptide import RandomGlycopeptideBuilder, forge_prediction_record

Sequence = sequence.Sequence
sequence_tokenizer_respect_sequons = sequence.sequence_tokenizer_respect_sequons
list_to_sequence = sequence.list_to_sequence
strip_modifications = sequence.strip_modifications
StubGlycopeptide = stub_glycopeptides.StubGlycopeptide


def n_unique_elements(seq):
    unique_elements = []
    for i in range(len(seq)):
        is_unique = True
        cur = seq[i]
        for pos in unique_elements:
            if compare_positions(cur, pos):
                is_unique = False
                break
        if is_unique:
            unique_elements.append(cur)
    return len(unique_elements)


def permute_sequence(seq):
    seq = copy.deepcopy(seq)
    ix_a = ix_b = 0
    length = len(seq) - 1
    total_perms = math.factorial(n_unique_elements(seq))
    perms_produced = 0
    while perms_produced < total_perms:
        ix_a = random.randint(0, length)
        ix_b = random.randint(0, length)
        if ix_a == ix_b:
            continue
        a = seq[ix_a]
        b = seq[ix_b]
        if compare_positions(a, b):
            continue
        seq[ix_a] = b
        seq[ix_b] = a
        perms_produced += 1
        yield copy.deepcopy(seq)

def compare_positions(a, b):
    # Are the residues are the same and the number of modifications are equal?
    res = (a[0] == b[0]) and (len(a[1]) == len(b[1]))
    if not res:
        return res
    # Compare each modification
    for i in range(len(a[1])):
        res = (a[1][i] == b[1][i])
        if not res:
            return res
    return res


def clean_tokenizer(seq):
    tokens, mods, glycan, n_term, c_term = sequence.sequence_tokenizer(seq)
    return tokens


def edit_distance(new_seq, prev_seq):
    previous = range(len(clean_tokenizer(prev_seq)) + 1)
    for i, new_pos in enumerate(clean_tokenizer(new_seq)):
        current = [i + 1]
        for j, prev_pos in enumerate(clean_tokenizer(prev_seq)):
            insertions = previous[j + 1] + 1
            deletions = current[j] + 1
            substitutions = previous[j] + (not compare_positions(new_pos, prev_pos))
            current.append(min(insertions, deletions, substitutions))
        previous = current
    return previous[-1]


def build_shuffle_sequences(ix_prediction, count=20, prefix_len=0, suffix_len=0,
                            iter_max=None, builder=None, random_only=False):
    '''
        :type scored_matches: PredictionResults
        :param scored_matches: pd.DataFrame having scored MS2 matches
        :param count: int > 0, number of random peptides to generate per match
        :param prefix_len: int >= 0, number of AA to preserve order on of the
                           start of the match sequence
        :param suffix_len: int >= 0, number of AA to preserve order on of the
                           end of the match sequence
    '''
    pname = multiprocessing.current_process().name
    iter_max = count * 10 if iter_max is None else iter_max
    shuffles = list()
    total_missing = 0
    ix, row = ix_prediction
    solutions = set()
    seq = Sequence(row.Glycopeptide_identifier).get_sequence(include_glycan=False)
    decoy_logger.info("[%s] Building %d decoys for  %s", pname, count, (ix, row.Glycopeptide_identifier, row.Calc_mass))
    solutions.add(seq)
    if not random_only:
        iter_count = 0
        clone = sequence_tokenizer_respect_sequons(row.Glycopeptide_identifier)
        pref = clone[:prefix_len]
        suf = clone[-suffix_len:]
        body = clone[prefix_len:(-suffix_len)]
        n_unique = n_unique_elements(body)
        min_diff = n_unique/3.0
        while(len(solutions) - 1 < count and iter_count < iter_max):
            # Permutations are very similar to the original sequence so first
            # transforming the sequence by shuffling it at random or reversing
            # it produce more heterogenous decoys
            # random.shuffle(body)
            body = body[::-1]
            prev = str(list_to_sequence(pref + list(body) + suf))
            solutions.add(prev)
            for shuffle in permute_sequence(body):  #itertools.permutations(body):
                clone = pref + list(shuffle) + suf
                res = str(list_to_sequence(clone))
                if edit_distance(res, prev) > (min_diff):
                    solutions.add(res)
                assert len(res) == len(seq)
                if(len(solutions) - 1) >= count:
                    break
            iter_count += 1
    solutions.discard(seq)
    decoys = []
    for shuffle in solutions:
        d = row.copy()
        d.Glycopeptide_identifier = shuffle + row.Glycan
        d._batch_id = ix
        decoys.append(d)

    short = count - len(solutions)
    if(short > 0):
        if not random_only:
            decoy_logger.info("%s was short %d sequences", str(seq), short)
        randomized = builder.generate_random(row.Calc_mass, short)
        if len(randomized) != short:
            decoy_logger.warning("(Glycopeptide: %s, Mass: %e) Was short %d decoys. Randomly generated only %d",
                                 row.Glycopeptide_identifier, row.Calc_mass, short, len(randomized))
            total_missing += short - len(randomized)
        for rand in randomized:
            forgery = forge_prediction_record(rand, row)
            forgery["_batch_id"] = ix
            decoys.append(forgery)

    decoy_search_space = []
    for decoy in decoys:
        try:
            frags = random_glycopeptide_to_fragments(decoy)
            frags["progenitor"] = row.Glycopeptide_identifier
            decoy_search_space.append(frags)
        except Exception, e:
            decoy_logger.error(e, exc_info=e)
            decoy_logger.warning("Decoy %s failed to be forged", decoy.Glycopeptide_identifier)
            raise
    decoy_logger.info("Finished forging %d decoys for %s", len(decoy_search_space), row.Glycopeptide_identifier)
    shuffles.append((decoy_search_space, total_missing))
    return shuffles


def random_glycopeptide_to_fragments(sequence_record):
    try:
        seq_obj = Sequence(sequence_record.Glycopeptide_identifier)
    except:
        print(sequence_record)
        raise
    glycan_map = {}
    modifications = []
    loc_rules = modification.get_position_modifier_rules_dict(seq_obj)
    for i, (aa, mods) in enumerate(seq_obj):
        for mod in mods:
            if "Glycan" in mod.name or "HexNAc" in mod.name:
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

    # Remove glycans from the sequence string to conform to the SequenceSpace expectations
    for site, glycan in glycan_map.items():
        # Don't discard anonymous HexNAcs. Downstream functions can handle them
        if glycan != "HexNAc":
            seq_obj.drop_modification(site, glycan)
    glycosylation_sites = glycan_map.keys()
    if sequence_record.Glycan is None:
        # Build the semicolon separated string for glycan compositions
        glycan_composition = []
        glycan_composition = [map(int, glycan.replace("Glycan", '').replace("[", "").replace("]", "").split(";"))
                              for glycan in glycan_map.values()]
        glycan_composition = map(sum, zip(*glycan_composition))
        glycan_composition_restring = "[" + ";".join(map(str, glycan_composition)) + "]"
    else:
        glycan_composition_restring = sequence_record.Glycan
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
    pep_stubs = StubGlycopeptide(peptide_seq, glycan_composition_restring,
                                 len(glycosylation_sites), glycan_composition_restring)
    stub_ions = pep_stubs.get_stubs()
    assert len(stub_ions) > 1
    oxonium_ions = pep_stubs.get_oxonium_ions()
    ions = {
        "MS1_Score": sequence_record.MS1_Score,
        "Obs_Mass": sequence_record.Obs_Mass,
        "Calc_mass": sequence_record.Calc_mass,
        "ppm_error": sequence_record.ppm_error,
        "Peptide": peptide_seq,
        "Peptide_mod": sequence_record.Peptide_mod,
        "Glycan": glycan_composition_restring,
        "vol": sequence_record.vol,
        "glyco_sites": len(glycan_map),
        "startAA": None,
        "endAA": None,
        "Seq_with_mod": seq_obj.get_sequence(include_glycan=False),
        "bare_b_ions": b_ions,
        "b_ions_with_HexNAc": b_ions_HexNAc,
        "bare_y_ions": y_ions,
        "y_ions_with_HexNAc": y_ions_HexNAc,
        "pep_stub_ions": stub_ions,
        "Oxonium_ions": oxonium_ions,
        "Glycopeptide_identifier": seq_obj.get_sequence(include_glycan=False) + glycan_composition_restring,
        "_batch_id": int(sequence_record._batch_id)
    }
    return ions


def taskmain(predictions_path, prefix_len=0, suffix_len=0,
             count=20, random_only=False, n_processes=5, out=None):

    start_time = time.time()
    predictions = prediction_tools.prepare_model_file(predictions_path)
    decoy_path = (predictions_path if out is None else out).rsplit("scored", 1)[0] + "decoy.db"
    metadata = predictions.metadata
    metadata = copy.deepcopy(metadata)
    metadata["tag"] = (metadata["tag"] if metadata["tag"] not in ["", None] else "") + "decoy"
    metadata["decoy_ratio"] = count

    decoy_logger.info("Saving metadata")
    metadata_store = sqlitedict.SqliteDict(decoy_path, tablename="metadata", flag='n')
    metadata_store.update(metadata)
    metadata_store.commit()
    enzyme = predictions.metadata.get("enzyme")
    starts = sum([enz["cleavage_start"] for enz in enzyme], [])
    ends = sum([enz["cleavage_end"] for enz in enzyme], [])
    decoy_logger.info("Building random sequence builder")
    decoy_logger.info("Parameters: {}".format(dict(prefix_len=prefix_len,
                                                   suffix_len=suffix_len, count=count,
                                                   random_only=random_only)))
    modification_table = modification.ModificationTable.bootstrap()
    builder = RandomGlycopeptideBuilder(
        ppm_error=predictions.metadata.get("ms1_ppm_tolerance", 10e-6),
        constant_modifications=map(modification_table.__getitem__,
                                   predictions.metadata.get("constant_modifications", [])),
        variable_modifications=map(modification_table.__getitem__,
                                   predictions.metadata.get("variable_modifciations", [])),
        glycans=glycans_from_predictions(predictions),
        cleavage_start=starts,
        cleavage_end=ends)
    task_fn = functools.partial(build_shuffle_sequences, prefix_len=prefix_len,
                                suffix_len=suffix_len, count=count,
                                builder=builder, random_only=random_only)
    decoy_logger.info("Beginning generation")
    decoy_sequence_store = sqlitedict.SqliteDict(decoy_path,
                                                 tablename="theoretical_search_space", flag='w')
    total_missing = 0
    preds_missed = []
    pool = None
    if n_processes > 1:
        pool = multiprocessing.Pool(n_processes)
        dispatcher = functools.partial(pool.imap_unordered, chunksize=1000)
        decoy_logger.info("Building decoys concurrently on %d processes", n_processes)
    else:
        decoy_logger.info("Building decoys sequentially")
        dispatcher = itertools.imap
    i = 0
    pred_iter = 0
    for decoy_tuple in (
        itertools.chain.from_iterable(
            dispatcher(task_fn, predictions.iterrows()))):
        decoys, missing = decoy_tuple
        for decoy in decoys:
            decoy_sequence_store[i] = decoy
            i += 1
            if i % 10000 == 0:
                decoy_sequence_store.commit()
        total_missing += missing
        if missing > 0:
            decoy_logger.info("%d missing", missing)
            preds_missed.append((predictions.ix[pred_iter].Glycopeptide_identifier, missing))
        pred_iter += 1
    if pool is not None:
        pool.close()
        pool.join()
    decoy_sequence_store.commit()
    decoy_sequence_store.close()
    metadata_store["total_missing_decoys"] = total_missing
    metadata_store["predictions_missing_decoys"] = preds_missed
    metadata_store.commit()
    metadata_store.close()
    decoy_logger.info("Done. %d Sequences missing. %f seconds elapsed",
                      total_missing, time.time() - start_time)
    return decoy_path


app = argparse.ArgumentParser("make_decoys")
app.add_argument("predictions_path", help="File path to the prediction results to make decoys for")
app.add_argument("-r", "--random-only", action="store_true")


if __name__ == '__main__':
    import sys
    logging.basicConfig(level="INFO")
    taskmain(sys.argv[1], 0, 1, 20, 1)
