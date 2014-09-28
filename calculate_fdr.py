import csv
import itertools

import random_glycopeptide
import classify_matches
import match_ions


def filter_fn(operating_data):
    ind = (operating_data.numStubs > 0) & \
        (operating_data.meanCoverage > 0.0) & \
        (operating_data.meanHexNAcCoverage > 0.0) & \
        (operating_data.MS2_Score > 0.5)
    return ind


def load_decoys_from_fasta(fasta_file_handle):
    '''Reads decoys from a FASTA format file, yielding them as an iterator'''
    defline = ''
    seq = ''
    try:
        while(True):
            # defline is unused
            defline = fasta_file_handle.next()
            seq = fasta_file_handle.next().replace("\n", "").replace("\r", "")
            yield dict(defline=defline, sequence=seq)
    except StopIteration:
        pass


def generate_theoretical_ions(fasta_file_path):
    ions = []
    for decoy in load_decoys_from_fasta(open(fasta_file_path)):
        ions.extend(random_glycopeptide.random_glycopeptide_to_sequence_space(decoy["sequence"]))
    return ions


def main(decoy_fasta, decon_data, ms1_tolerance=1e-5, ms2_tolerance=2e-5,
         scored_matches_count=0, num_decoys_per_real_mass=20.0):
        tempfile_name = "___temp_decoy_ion_space.csv"
        decoy_ions = generate_theoretical_ions(decoy_fasta)
        fh = open("___temp_decoy_ion_space.csv", "wb")
        w = csv.DictWriter(fh, decoy_ions[0].keys())
        w.writeheader()
        w.writerows(decoy_ions)
        match_file = match_ions.match_frags(tempfile_name, decon_data, ms1_tolerance, ms2_tolerance)
        matches = [m for m in csv.DictReader(open(match_file))]
        num_decoy_matches = len(matches)
        print(num_decoy_matches)
        total_assignments = scored_matches_count + num_decoy_matches
        fdr = (num_decoy_matches / float(total_assignments)) * (1 + (1 / num_decoys_per_real_mass))
        return fdr


if __name__ == '__main__':
    import argparse
    app = argparse.ArgumentParser()
    app.add_argument("decoy_fasta")
    app.add_argument("decon_data")
    app.add_argument("-s", "--scored-matches", required=True,
                     help="Path to results of matching or a number indicating the number of 'real' matches")

    args = app.parse_args()

    args = args.__dict__
    scored = args.pop("scored_matches")
    try:
        count = int(scored)
    except:
        try:
            import pandas as pd
            scored_table = pd.read_csv(scored)
            count = len(scored_table.index)
        except:
            print("Could not determine the number of matches from %s.\
                Could not interpret as a number or file path" % scored)
            exit(-1)
    args['scored_matches_count'] = count
    main(**args)
