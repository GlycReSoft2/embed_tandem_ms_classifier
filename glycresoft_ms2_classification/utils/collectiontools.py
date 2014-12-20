from collections import defaultdict
import itertools


def groupby(ungrouped_list, key_fn=lambda x: x, kind=list):
    if kind == list:
        return groupby_list(ungrouped_list, key_fn)
    elif kind == set:
        return groupby_set(ungrouped_list, key_fn)
    else:
        raise TypeError("Collection {kind} is not defined".format(kind=kind))


def groupby_list(ungrouped_list, key_fn):
    groups = defaultdict(list)
    for item in ungrouped_list:
        key_value = key_fn(item)
        groups[key_value].append(item)
    return groups


def groupby_set(ungrouped_list, key_fn):
    groups = defaultdict(set)
    for item in ungrouped_list:
        key_value = key_fn(item)
        groups[key_value].add(item)
    return groups


def _repeat_label(group_size):
    label = 0
    while(True):
        for value in itertools.repeat(label, group_size):
            yield value
        label += 1


def chunk_iterator(iterator, chunk_size):
    labeler = _repeat_label(chunk_size)
    for label, chunk in itertools.groupby(iterator, lambda x: labeler.next()):
        yield list(chunk)
