import multiprocessing
import functools
import itertools
import json
import copy

from glycresoft_ms2_classification import prediction_tools
from glycresoft_ms2_classification import calculate_fdr

shuffle_sequence = calculate_fdr.shuffle_sequence
random_glycopeptide_to_fragments = calculate_fdr.random_glycopeptide_to_fragments


def repackage_shuffle(prediction, shuffles):
    decoys = [None] * len(shuffles)
    for i, shuffle in enumerate(shuffles):
        d = prediction.copy()
        d.Glycopeptide_identifier = shuffle + prediction.Glycan
        decoys[i] = d
    return decoys


def make_decoys(ix_prediction, prefix_len, suffix_len, count):
    ix, prediction = ix_prediction
    shuffles = shuffle_sequence(prediction.Glycopeptide_identifier, prefix_len, suffix_len, count, count * 10)
    decoys = repackage_shuffle(prediction, shuffles)
    for decoy in decoys:
        decoy["_batch_id"] = ix
    decoy_dicts = map(random_glycopeptide_to_fragments, decoys)
    return decoy_dicts


def main(predictions_path, prefix_len, suffix_len, count, n_processes):
    predictions = prediction_tools.prepare_model_file(predictions_path)
    metadata = predictions.metadata
    metadata = copy.deepcopy(metadata)
    metadata["tag"] = (metadata["tag"] if metadata["tag"] not in {"", None} else "") + "decoy"
    metadata["decoy_ratio"] = count
    pool = multiprocessing.Pool(n_processes)
    task_fn = functools.partial(make_decoys, prefix_len=prefix_len, suffix_len=suffix_len, count=count)
    decoy_fragments = {
        "metadata": metadata,
        "theoretical_search_space": []
    }
    acc = []
    for decoy in itertools.chain.from_iterable(pool.imap_unordered(task_fn, predictions.iterrows(), 100)):
        acc.append(decoy)
    decoy_fragments["theoretical_search_space"] = acc
    print("Serializing")
    decoy_file = open(predictions_path[:-4] + "decoy.ion_space.json", "wb")

    for chunk in json.JSONEncoder().iterencode(decoy_fragments):
        decoy_file.write(chunk)

    decoy_file.close()
    print(decoy_file.name)

if __name__ == '__main__':
    import sys
    main(sys.argv[1], 0, 1, 20, 6)
