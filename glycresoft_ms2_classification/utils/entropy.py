import math
from collections import Counter

from structure.sequence import Sequence

def count_positions(sequence):
    counts = Counter()
    for pos in sequence:
        ipos = (pos[0].name, tuple(map(lambda x: x.name, pos[1])))
        counts[ipos] += 1
    return counts

def relative_frequency(counts):
    freqs = dict()
    total = float(sum(counts.values()))
    for ipos, count in counts.items():
        freqs[ipos] = count/total
    return freqs

def sequence_entropy(sequence):
    counts = count_positions(sequence)
    freqs = relative_frequency(counts)
    entropy = 0.0
    for ipos, freq in freqs.items():
        entropy += freq * math.log(1/freq, 2)
    return entropy

def main(sequences):
    sequences = map(lambda x: Sequence(x.split("[")[0]), sequences)
    entropies = map(sequence_entropy, sequences)
    return entropies

if __name__ == '__main__':
    import sys
    import os
    sequence_file = sys.argv[1]
    seqs = open(sequence_file).readlines()
    entropies = main(seqs)
    from matplotlib import pyplot as plt
    from pandas import Series
    es = Series(entropies)
    es.hist()
    plt.savefig(os.path.splitext(sequence_file)[0] + "_entropy_hist.png")
