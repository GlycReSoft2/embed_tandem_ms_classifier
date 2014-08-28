
#! /usr/bin/python

# Some of the codes are modified based on KK's scripts.

# Input:
# Arg 1 -- csv file containing glycopeptide composition
# Arg 2 -- list specifying the candidate glycosylation sites.
# Example: python HH_glycopeptide.py test_data.csv sample_glycosites.txt >
# report.txt

import csv
#import json
import os
import re
#import sys

from collections import defaultdict

#from structure.residue import Residue
from structure.modification import RestrictedModificationTable
from structure.modification import ModificationTable
from structure.sequence_space import SequenceSpace
#from structure.sequence import Sequence
from structure.stub_glycopeptides import stubs
#from structure.fragment import Fragment

KEYS = [
    "MS1_Score", "Obs_Mass", "Calc_mass", "ppm_error", "Peptide", "Peptide_mod", "Glycan", "vol",
    "glyco_sites", "startAA", "endAA", "Seq_with_mod", "Glycopeptide_identifier", "Oxonium_ions",
    "pep_stub_ions", "bare_b_ions", "bare_y_ions", "b_ions_with_HexNAc", "y_ions_with_HexNAc"]


class TheoreticalIonFragment(object):

    mod_pattern = re.compile(r'(\d+)(\w+)')

    def __init__(self, **fields):
        data_dict = defaultdict(lambda: None)
        data_dict.update(fields)
        self.ms1_score = data_dict[KEYS[0]]
        self.obs_mass = data_dict[KEYS[1]]
        self.calc_mass = data_dict[KEYS[2]]
        self.ppm_error = data_dict[KEYS[3]]
        self.peptide = data_dict[KEYS[4]]
        self.peptide_mod = data_dict[KEYS[5]]
        self.glycan = data_dict[KEYS[6]]
        self.vol = data_dict[KEYS[7]]
        self.glyco_sites = data_dict[KEYS[8]]
        self.start_aa = data_dict[KEYS[9]]
        self.end_aa = data_dict[KEYS[10]]
        self.seq_with_mod = data_dict[KEYS[11]]
        self.glycopeptide_identifier = data_dict[KEYS[12]]
        self.oxonium_ions = data_dict[KEYS[13]]
        self.pep_stub_ions = data_dict[KEYS[14]]
        self.bare_b_ions = data_dict[KEYS[15]]
        self.bare_y_ions = data_dict[KEYS[16]]
        self.b_ions_with_hexnac = data_dict[KEYS[17]]
        self.y_ions_with_hexnac = data_dict[KEYS[18]]

    @classmethod
    def get_peptide_modifications(cls, data, modification_table):
        items = cls.mod_pattern.findall(data)
        mod_list = []
        for i in items:
            mod = modification_table.get_modification(i[1], -1, int(i[0]))
            mod_list.append(mod)
        return mod_list

    @classmethod
    def get_search_space(cls, row, glycan_identities, glycan_sites, seq_str, mod_list):
        glycan_compo = {}
        for g in glycan_identities:
            glycan_compo[g] = int(row[''.join(['G:', g])])
        seq_space = SequenceSpace(
            seq_str, glycan_compo, glycan_sites, mod_list)
        # print(seq_space)
        return seq_space


def main(result_file, site_file, constant_modification_list=None, variable_modification_list=None, output_file=None):
    if output_file is None:
        output_file = os.path.splitext(result_file)[0] + '.theoretical_ions'

    modification_table = RestrictedModificationTable.bootstrap(
        constant_modification_list, variable_modification_list)

    if constant_modification_list is None and variable_modification_list is None:
        modification_table = ModificationTable.bootstrap()

    #print("Reading %s" % result_file)
    compo_dict = csv.DictReader(open(result_file, "r"), delimiter=",")
    site_list = [line.strip() for line in open(site_file, "r")]
    site_list = list(map(int, site_list))

    colnames = compo_dict.fieldnames
    match_pattern = re.compile(r'^G:(\w+)')
    glycan_identity = []
    for col in colnames:
        if match_pattern.match(col):
            glycan_identity.extend(match_pattern.findall(col))

    fragment_info = []
    # Each row as a dict.
    #print("Creating set of all theoretical glycopeptides")
    for i, row in enumerate(compo_dict):
        # if (i % 1000) == 0:
            #print("%d total" % i)
        # For each row, read the information of glycan compound, peptide sequence,
        # ofGlycanAttachmentToPeptide, PeptideModification
        score = row['Score']
        precur_mass = row['MassSpec MW']
        theo_mass = row['Hypothesis MW']
        glycan_comp = row['Compound Key']
        seq_str = row['PeptideSequence']
        pep_mod = row['PeptideModification']
        pep_mis = row['PeptideMissedCleavage#']
        num_sites = int(row['#ofGlycanAttachmentToPeptide'])
        mass_error = row['PPM Error']
        volume = row['Total Volume']

        if (seq_str == '') or (num_sites == 0):  # Not interested
            continue

        # Compute the set of modifications that can occur.
        mod_list = TheoreticalIonFragment.get_peptide_modifications(
            row['PeptideModification'], modification_table)

        # Get the start and end positions of fragment relative to the
        start_pos = int(row['StartAA'])
        end_pos = int(row['EndAA'])
        glycan_sites = set(site_list).intersection(
            range(start_pos, end_pos + 1))

        # No recorded sites, skip this component.
        if len(glycan_sites) == 0:
            continue

        # Adjust the glycan_sites to relative position
        glycan_sites = [x - start_pos for x in glycan_sites]

        ss = TheoreticalIonFragment.get_search_space(
            row, glycan_identity, glycan_sites, seq_str, mod_list)
        seq_list = ss.get_theoretical_sequence(num_sites)
        for s in seq_list:
            #print(seq_str, '\t', row['Compound Key'], mod_str, '\n', s.getSequence())
            seq_mod = s.getSequence()

            b_type = s.getFragments('B')
            b_ions = []
            b_ions_HexNAc = []
            for b in b_type:
                for fm in b:
                    key = fm.getFragmentName()
                    if key == "B1" or re.search(r'B1\+', key):
                        # B1 Ions aren't actually seen in reality, but are an artefact of the generation process
                        # so do not include them in the output
                        continue
                    mass = fm.getMass()
                    if "HexNAc" in key:
                        b_ions_HexNAc.append({"key": key, "mass": mass})
                    else:
                        b_ions.append({"key": key, "mass": mass})

            y_type = s.getFragments('Y')
            y_ions = []
            y_ions_HexNAc = []
            for y in y_type:
                for fm in y:
                    key = fm.getFragmentName()
                    mass = fm.getMass()
                    if "HexNAc" in key:
                        y_ions_HexNAc.append({"key": key, "mass": mass})
                    else:
                        y_ions.append({"key": key, "mass": mass})

            pep_stubs = stubs(seq_str, pep_mod, num_sites, glycan_comp)
            stub_ions = pep_stubs.getStubs()
            oxonium_ions = pep_stubs.getOxoniums()

            fragment_info.append(
                {"MS1_Score": score, "Obs_Mass": precur_mass, "Calc_mass": theo_mass, "ppm_error": mass_error,
                 "Peptide": seq_str, "Peptide_mod": pep_mod, "Glycan": glycan_comp, "vol": volume, "glyco_sites": num_sites,
                 "startAA": start_pos, "endAA": end_pos, "Seq_with_mod": seq_mod, "Glycopeptide_identifier": seq_mod + glycan_comp,
                 "Oxonium_ions": oxonium_ions, "pep_stub_ions": stub_ions, "bare_b_ions": b_ions, "bare_y_ions": y_ions,
                 "b_ions_with_HexNAc": b_ions_HexNAc, "y_ions_with_HexNAc": y_ions_HexNAc})

    fh = open(output_file + '.csv', 'wb')
    dict_writer = csv.DictWriter(fh, KEYS)
    dict_writer.writer.writerow(KEYS)
    dict_writer.writerows(fragment_info)
    fh.close()
    return output_file + '.csv'


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description="Generate all theoretical ions from input peptides and possible glycosylation sites")
    parser.add_argument("-f", "--fragment-file", type=str,
                        help="The csv file produced by Glycresoft, describing a peptide")
    parser.add_argument("-s", "--site-file", type=str,
                        help="A file listing each position along the peptide where ")
    parser.add_argument("-o", "--output-file", type=str, default=None,
                        help="The name of the file to output results to. Defaults to `fragment-file`_theoretical_ions")
    parser.add_argument(
        "--constant-modification-list", type=str, action="append", default=None,
        help="Pass the list of constant modifications to include in the sequence search space")
    parser.add_argument(
        "--variable-modification-list", type=str, action="append", default=None,
        help="Pass the list of variable modifications to include in the sequence search space")
    args = parser.parse_args()

    main(result_file=args.fragment_file, site_file=args.site_file,
         constant_modification_list=args.constant_modification_list,
         variable_modification_list=args.variable_modification_list,
         output_file=args.output_file)
