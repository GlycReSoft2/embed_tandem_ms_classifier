import csv
import re
import argparse
import sys
import os
import itertools

from pyteomics.parser import expasy_rules
from glycresoft_ms2_classification.structure import sequence, modification

Modification = modification.Modification
Sequence = sequence.Sequence


def try_cleave(row, rule, min_length=8, missed_cleavages=1):
    results = sequence.cleave(row["Sequence"], rule, min_length=min_length,
                              missed_cleavages=missed_cleavages)
    original = Sequence(row["Sequence"])
    original_glycan_sites = original.n_glycan_sequon_sites
    if len(results) == 1:
        yield row
    else:
        mods = None
        mod_counts = []
        if row["Modifications"].strip() != "":
            mod_list = [s for s in row["Modifications"].split(" ") if s != ""]
            mods = [None] * len(mod_list)
            mod_counts = [0] * len(mod_list)
            for i, mod in enumerate(mod_list):
                mod_count, mod_type = [t for t in re.split(r"(\d*)(.+)", mod) if t != ""]
                mod_count = int(mod_count) if len(mod_count) > 0 else 1
                mod = Modification(mod_type)
                mods[i] = mod
                mod_counts[i] = mod_count
        start = int(row["Start"])
        end = int(row["End"])
        glycan_sites_left = 0
        mods_left = 0
        for cleaved_seq, cleave_start, cleave_end in results:
            seq = Sequence(cleaved_seq)
            glycan_sites_left += len(seq.n_glycan_sequon_sites)
            mod_str = ""
            mass = seq.mass
            if mods is not None:
                for i, mod in enumerate(mods):
                    new_mod_count = len(mod.find_valid_sites(seq))
                    mods_left += new_mod_count
                    mod_str += str(new_mod_count)+mod_type if new_mod_count > 0 else ""
                    mass += (mod.mass * min(new_mod_count, mod_count))

            new_dict = dict(
                Select=row["Select"],
                Start=(start + cleave_start + 1),
                End=(end + cleave_end + 1),
                Sequence=cleaved_seq,
                Modifications=mod_str,
                Mass=mass)
            yield new_dict

        if glycan_sites_left < original_glycan_sites:
            yield row

app = argparse.ArgumentParser()
app.add_argument("gpe", type=str, nargs="?", default='--')
app.add_argument("-o", "--out", type=argparse.FileType('w'), default=None)
app.add_argument("enzyme_name", choices=expasy_rules.keys())
app.add_argument("-m", "--missed_cleavages", type=int, default=1)


def main():
    args = app.parse_args()
    gpe = args.gpe
    if gpe == "--":
        gpe = sys.stdin
    else:
        gpe = open(gpe)
    enzyme_name = args.enzyme_name
    missed_cleavages = args.missed_cleavages
    out = args.out if args.out is not None else sys.stdout

    rule = expasy_rules[enzyme_name]
    reader = csv.DictReader(gpe)
    writer = csv.DictWriter(out, fieldnames=reader.fieldnames)
    writer.writeheader()
    writer.writerows(
        itertools.chain.from_iterable(
            try_cleave(row, rule, missed_cleavages=missed_cleavages) for row in reader))
    sys.stdout.flush()

if __name__ == '__main__':
    main()
