import csv
import os
import re
import json
import multiprocessing
import logging
import functools
import itertools

from structure.modification import RestrictedModificationTable
from structure.modification import ModificationTable
from structure.sequence_space import SequenceSpace
from structure.stub_glycopeptides import StubGlycopeptide
from structure import constants

from proteomics import get_enzyme

mod_pattern = re.compile(r'(\d+)(\w+)')


def get_peptide_modifications(data, modification_table):
    items = mod_pattern.findall(data)
    mod_list = []
    for i in items:
        if i[1] == '':
            continue
        mod = modification_table.get_modification(i[1], -1, int(i[0]))
        mod_list.append(mod)
    return mod_list


def get_search_space(row, glycan_identities, glycan_sites, seq_str, mod_list):
    glycan_compo = {}
    for g in glycan_identities:
        glycan_compo[g] = int(row[''.join(['G:', g])])
    seq_space = SequenceSpace(
        seq_str, glycan_compo, glycan_sites, mod_list)
    return seq_space


def get_glycan_identities(colnames):
    glycan_identity = []
    extract_state = False
    for col in colnames:
        if col == "Hypothesis MW":
            extract_state = True
            continue
        elif col == "Adduct/Replacement":
            extract_state = False
            break
        elif extract_state:
            glycan_identity.append(col.replace("G:", ""))
    return glycan_identity


def generate_fragments(
    seq, num_sites=0, glycan_comp=None, pep_mod=None, seq_str="", theo_mass=0., mass_error=0.,
        score=0, precur_mass=0, volume=0, start_pos=-1, end_pos=-1):
    seq_mod = seq.get_sequence()
    b_type = seq.get_fragments('B')
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
            if "HexNAc" in key:
                b_ions_HexNAc.append({"key": key, "mass": mass})
            else:
                b_ions.append({"key": key, "mass": mass})

    y_type = seq.get_fragments('Y')
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

    pep_stubs = StubGlycopeptide(seq_str, pep_mod, num_sites, glycan_comp)
    stub_ions = pep_stubs.get_stubs()
    oxonium_ions = pep_stubs.get_oxonium_ions()

    return {
        "MS1_Score": score, "Obs_Mass": precur_mass, "Calc_mass": theo_mass, "ppm_error": mass_error,
        "Peptide": seq_str, "Peptide_mod": pep_mod, "Glycan": glycan_comp, "vol": volume, "glyco_sites": num_sites,
        "startAA": start_pos, "endAA": end_pos, "Seq_with_mod": seq_mod,
        "Glycopeptide_identifier": seq_mod + glycan_comp,
        "Oxonium_ions": oxonium_ions, "pep_stub_ions": stub_ions, "bare_b_ions": b_ions, "bare_y_ions": y_ions,
        "b_ions_with_HexNAc": b_ions_HexNAc, "y_ions_with_HexNAc": y_ions_HexNAc}


def process_predicted_ms1_ion(row, modification_table, site_list, glycan_identity):
    score = float(row['Score'])
    precur_mass = float(row['MassSpec MW'])
    theo_mass = float(row['Hypothesis MW'])
    glycan_comp = row['Compound Key']
    seq_str = row['PeptideSequence']
    pep_mod = row['PeptideModification']
    # pep_mis = row['PeptideMissedCleavage#']
    num_sites = int(row['#ofGlycanAttachmentToPeptide'])
    mass_error = float(row['PPM Error'])
    volume = float(row['Total Volume'])

    if (seq_str == '') or (num_sites == 0):
        return []

    # Compute the set of modifications that can occur.
    mod_list = get_peptide_modifications(
        row['PeptideModification'], modification_table)

    # Get the start and end positions of fragment relative to the
    start_pos = int(row['StartAA'])
    end_pos = int(row['EndAA'])
    glycan_sites = set(site_list).intersection(
        range(start_pos, end_pos + 1))

    # No recorded sites, skip this component.
    if len(glycan_sites) == 0:
        return []

    # Adjust the glycan_sites to relative position
    glycan_sites = [x - start_pos for x in glycan_sites]
    ss = get_search_space(
        row, glycan_identity, glycan_sites, seq_str, mod_list)
    seq_list = ss.get_theoretical_sequence(num_sites)
    fragments = [generate_fragments(seq, num_sites, glycan_comp, pep_mod, seq_str,
                                    theo_mass, mass_error, score, precur_mass,
                                    volume, start_pos, end_pos)
                 for seq in seq_list]
    return fragments


def main(result_file, site_file, constant_modification_list=None, variable_modification_list=None,
         enzyme_info=None, n_processes=4, output_file=None):
    if output_file is None:
        output_file = os.path.splitext(result_file)[0] + '.theoretical_ions'
    modification_table = RestrictedModificationTable.bootstrap(
        constant_modification_list, variable_modification_list)

    if constant_modification_list is None and variable_modification_list is None:
        modification_table = ModificationTable.bootstrap()

    if isinstance(site_file, basestring):
        site_list = [line.strip() for line in open(site_file, "r")]
        site_list = list(map(int, site_list))
    else:
        site_list = site_file

    compo_dict = csv.DictReader(open(result_file, "r"), delimiter=",")
    colnames = compo_dict.fieldnames
    glycan_identity = get_glycan_identities(colnames)
    enzyme_info = map(get_enzyme, enzyme_info)

    metadata = {
        "glycan_identities": glycan_identity,
        "constant_modifications": constant_modification_list,
        "variable_modifications": variable_modification_list,
        "site_list": site_list,
        "ms1_output_file": result_file,
        "enzyme": enzyme_info,
        "enable_partial_hexnac_match": constants.PARTIAL_HEXNAC_LOSS
    }

    print(constants.PARTIAL_HEXNAC_LOSS)

    fragment_info = []
    pool = multiprocessing.Pool(n_processes)

    task_fn = functools.partial(process_predicted_ms1_ion, modification_table=modification_table,
                                site_list=site_list, glycan_identity=glycan_identity)

    fragment_info = list(
        itertools.chain.from_iterable(pool.imap(task_fn, compo_dict, chunksize=1000)))

    if output_file is not False:
        print("Writing to file")
        fh = open(output_file + '.json', 'wb')
        data = {
            "metadata": metadata,
            "theoretical_search_space": fragment_info
        }
        json.dump(data, fh)
        fh.close()
    pool.close()
    pool.join()
    return output_file + '.json', data


def taskmain():
    import argparse
    parser = argparse.ArgumentParser(
        description="Generate all theoretical ions from input peptides and possible glycosylation sites")
    parser.add_argument("-f", "--fragment-file", type=str,
                        help="The csv file produced by Glycresoft, describing a peptide")
    parser.add_argument("-s", "--site-file", type=str,
                        help="A file listing each position along the peptide where glycosylation is predicted")
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
