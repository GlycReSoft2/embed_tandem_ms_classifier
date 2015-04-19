import csv
import os
import re
import datetime
import multiprocessing
import logging

import functools
import itertools

from sqlitedict import SqliteDict

from structure.modification import RestrictedModificationTable
from structure.modification import ModificationTable
from structure.sequence_space import SequenceSpace
from structure.stub_glycopeptides import StubGlycopeptide
from structure import constants

from proteomics import get_enzyme

logger = logging.getLogger("search_space_builder")
mod_pattern = re.compile(r'(\d+)(\w+)')
g_colon_prefix = "G:"


def _get_glycan_counts(mapping, glycan_identities):
    """Given a mapping and list of glycan names, get the mapped integer counts
    for each glycan
    """
    try:
        return {g: int(mapping[g_colon_prefix + g]) for g in glycan_identities}
    except:
        global g_colon_prefix
        g_colon_prefix = ""
        return {g: int(mapping[g_colon_prefix + g]) for g in glycan_identities}


class MS1GlycopeptideResult(object):
    '''
    Describes a row ofthe MS1 Results with format-agnostic basic
    construction and a CSV specific mapping
    '''
    @classmethod
    def from_csvdict(cls, glycan_identities=None, **kwargs):
        if glycan_identities is None:
            glycan_identities = ()
        score = float(kwargs.get("Score", 0.))
        theoretical_precursor_mass = float(kwargs.get("Hypothesis MW"))
        observed_precursor_mass = float(kwargs.get("MassSpec MW"))
        glycan_composition_str = kwargs.get("Compound Key")
        peptide_sequence = kwargs.get("PeptideSequence")
        peptide_modifications = kwargs.get("PeptideModification")
        num_missed_cleavages = int(kwargs.get("PeptideMissedCleavage"))
        num_glycosylation_sites = int(kwargs.get("#ofGlycanAttachmentToPeptide"))
        mass_error = float(kwargs.get("PPM Error"))
        volume = float(kwargs.get("Total Volume"))
        start_pos = int(kwargs.get("StartAA"))
        end_pos = int(kwargs.get("EndAA"))
        glycan_composition_map = _get_glycan_counts(kwargs, glycan_identities)
        return cls(
            score=score, theoretical_precursor_mass=theoretical_precursor_mass,
            observed_precursor_mass=observed_precursor_mass, glycan_composition_str=glycan_composition_str,
            peptide_sequence=peptide_sequence, peptide_modifications=peptide_modifications,
            num_missed_cleavages=num_missed_cleavages, num_glycosylation_sites=num_glycosylation_sites,
            mass_error=mass_error, volume=volume, start_pos=start_pos, end_pos=end_pos,
            glycan_composition_map=glycan_composition_map)

    def __init__(self, score=None, theoretical_precursor_mass=None,
                 observed_precursor_mass=None, glycan_composition_str=None,
                 peptide_sequence=None, peptide_modifications=None, num_missed_cleavages=None,
                 num_glycosylation_sites=None, mass_error=None, volume=None, start_pos=None,
                 end_pos=None, glycan_composition_map=None):
        self.score = score
        self.theoretical_precursor_mass = theoretical_precursor_mass
        self.observed_precursor_mass = observed_precursor_mass
        self.glycan_composition_str = glycan_composition_str
        self.peptide_sequence = peptide_sequence
        self.peptide_modifications = peptide_modifications
        self.num_missed_cleavages = num_missed_cleavages
        self.num_glycosylation_sites = num_glycosylation_sites
        self.mass_error = mass_error
        self.volume = volume
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.glycan_composition_map = glycan_composition_map

    def as_dict(self):
        dict_form = {
            "MS1_Score": self.score,
            "Obs_Mass": self.observed_precursor_mass,
            "Calc_mass": self.theoretical_precursor_mass, "ppm_error": self.mass_error,
            "Peptide": self.peptide_sequence, "Peptide_mod": self.peptide_modifications,
            "Glycan": self.glycan_composition_str, "vol": self.volume, "glyco_sites": self.num_sites,
            "startAA": self.start_pos, "endAA": self.end_pos
        }
        return dict_form


def get_glycan_identities(csv_columns):
    """Extract glycan names from the column headers of MS1 Results CSV files

    Parameters
    ----------
    csv_columns: list
        A list of strings containing column names from the MS1 Results.

    Returns
    -------
    glycan_identities : list
        The list of glycan identities
    """
    glycan_identity = []
    extract_state = False
    for col in csv_columns:
        if col == "Hypothesis MW":
            extract_state = True
            continue
        elif col == "Adduct/Replacement":
            extract_state = False
            break
        elif extract_state:
            glycan_identity.append(col.replace("G:", ""))
    logger.info("Glycan identities found: %s", glycan_identity)
    return glycan_identity


def get_peptide_modifications(peptide_mod_str, modification_table):
    """Maps the Protein Prospector modifications to :class:`.modifications.Modification`
    objects through `modification_table`

    Parameters
    ----------
    peptide_mod_str : str
        A string containing peptide modifications from ProteinProspector in the MS1 CSV
    modification_table : :class:`.modifications.ModificationTable`
        A mapping of possible :class:`.modifications.ModificationRule` functions provided
        by the user.

    Returns
    -------
    :class:`list` of :class:`.modifications.Modification` objects
    """
    items = mod_pattern.findall(peptide_mod_str)
    mod_list = []
    for i in items:
        if i[1] == '':
            continue
        mod = modification_table.get_modification(i[1], -1, int(i[0]))
        mod_list.append(mod)
    return mod_list


def get_search_space(ms1_result, glycan_sites, mod_list):
    """Create a :class:`.sequence_space.SequenceSpace` object from the collection
    of information interpolated from the MS1 Results

    Parameters
    ----------
    ms1_result : :class:`MS1GlycopeptideResult`
        The row of MS1 Results being operated on
    glycan_sites : list
        A list of sites along the parent sequence that are glycosylated
    mod_list : list
        The list of present modifications

    Returns
    -------
    :class:`.sequence_space.SequenceSpace`
        The generator of possible :class:`.sequence.Sequence` objects
    """

    seq_space = SequenceSpace(
        ms1_result.peptide_sequence,
        ms1_result.glycan_composition_map,
        glycan_sites, mod_list)
    return seq_space


def generate_fragments(seq, ms1_result):
    """Consumes a :class:`.sequence.Sequence` object, and the contents of an MS1 Result row to
    generate the set of all theoretically observed fragments

    Parameters
    ----------
    seq: :class:`.sequence.Sequence`
        The binding of modifications to particular sites on a peptide sequence
    ms1_result: :class:`MS1GlycopeptideResult`
        Description of the precursor match

    Returns
    -------
    dict:
        Collection of theoretical ions from the given sequence,
        as well as the precursor information.
    """
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

    pep_stubs = StubGlycopeptide(
        ms1_result.peptide_sequence,
        ms1_result.peptide_modifications,
        ms1_result.num_sites,
        ms1_result.glycan_composition_str)

    stub_ions = pep_stubs.get_stubs()
    oxonium_ions = pep_stubs.get_oxonium_ions()

    fragments_dict = {
        "Seq_with_mod": seq_mod,
        "Glycopeptide_identifier": seq_mod + ms1_result.glycan_composition_str,
        "Oxonium_ions": oxonium_ions,
        "pep_stub_ions": stub_ions,
        "bare_b_ions": b_ions,
        "bare_y_ions": y_ions,
        "b_ions_with_HexNAc": b_ions_HexNAc,
        "y_ions_with_HexNAc": y_ions_HexNAc
    }

    fragments_dict.update(ms1_result.as_dict())

    return fragments_dict


def process_predicted_ms1_ion(row, modification_table, site_list, glycan_identity):
    """Multiprocessing dispatch function to generate all theoretical sequences and their
    respective fragments from a given MS1 result

    Parameters
    ----------
    row: dict
        Line mapping from an MS1 results csv file
    modification_table: :class:`.modifications.RestrictedModificationTable`
        Modification table limited to only the rules specified by the user
    site_list: list of int
        List of putative glycosylation sites from the parent protein
    glycan_identity: list
        List of glycan or monosaccaride names

    Returns
    -------
    list of dicts:
        List of each theoretical sequence and its fragment ions
    """
    ms1_result = MS1GlycopeptideResult.from_csvdict(glycan_identity, **row)

    if (ms1_result.peptide_sequence == '') or (ms1_result.num_glycosylation_sites == 0):
        return []

    # Compute the set of modifications that can occur.
    mod_list = get_peptide_modifications(
        ms1_result.peptide_modifications, modification_table)

    # Get the start and end positions of fragment relative to the
    glycan_sites = set(site_list).intersection(
        range(ms1_result.start_pos, ms1_result.end_pos + 1))

    # No recorded sites, skip this component.
    if len(glycan_sites) == 0:
        return []

    # Adjust the glycan_sites to relative position
    glycan_sites = [x - ms1_result.start_pos for x in glycan_sites]
    ss = get_search_space(
        ms1_result, glycan_sites, mod_list)
    seq_list = ss.get_theoretical_sequence(ms1_result.num_sites)
    fragments = [generate_fragments(seq, ms1_result)
                 for seq in seq_list]
    return fragments


class TheoreticalSearchSpace(object):
    '''
    Describe the process of generating all theoretical sequences and their fragments
    from an MS1 Results CSV, a collection of constant and variable peptide modifications,
    and a list of putative glycosylation_sites.

    Uses an :class:`sqlitedict.SqliteDict` instance as a storage backend.

    Includes a single- and multi-process compatible implementation. The more processes used,
    the more memory must be allocated to buffer results.
    '''
    table_name = "theoretical_search_space"

    def __init__(self, ms1_results_file, db_file_name=None,
                 glycosylation_sites=None, enzyme_info=None,
                 site_list=None,
                 constant_modifications=None,
                 variable_modifications=None, n_processes=4):
        if db_file_name is None:
            db_file_name = os.splitext(ms1_results_file)[0] + '.db'
        self.store = SqliteDict(db_file_name, tablename=self.table_name, journal_mode="OFF")
        self.metadata = SqliteDict(db_file_name, tablename='metadata', journal_mode="OFF")
        self.ms1_results_file = ms1_results_file

        modification_table = RestrictedModificationTable.bootstrap(constant_modifications, variable_modifications)
        if constant_modifications is None and variable_modifications is None:
            modification_table = ModificationTable.bootstrap()
        self.modification_table = modification_table
        self.glycosylation_site = site_list
        self.ms1_results_reader = csv.DictReader(open(self.ms1_results_file))

        self.n_processes = n_processes

        self.glycan_identities = get_glycan_identities(self.ms1_results_reader.colnames)
        enzyme_info = map(get_enzyme, enzyme_info)
        tag = datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d-%H%M%S")
        self.metadata.update({
            "glycan_identities": self.glycan_identities,
            "enzyme_info": enzyme_info,
            "site_list": site_list,
            "constant_modification_list": constant_modifications,
            "variable_modification_list": variable_modifications,
            "ms1_output_file": ms1_results_file,
            "enzyme": enzyme_info,
            "tag": tag,
            "enable_partial_hexnac_match": constants.PARTIAL_HEXNAC_LOSS
        })
        self.metadata.commit()
        self.metadata.close()

    def run(self):
        '''
        Execute the algorithm on :attr:`n_processes` processes
        '''
        task_fn = functools.partial(process_predicted_ms1_ion, modification_table=self.modification_table,
                                    site_list=self.site_list, glycan_identity=self.glycan_identities)
        results_store = self.store
        cntr = 0
        if self.n_processes > 1:
            worker_pool = multiprocessing.Pool(self.n_processes)
            logger.debug("Building theoretical sequences concurrently")
            for res in (itertools.chain.from_iterable(
                    worker_pool.imap(task_fn, self.ms1_results_reader, chunksize=500))):
                results_store[cntr] = res
                cntr += 1
                if (cntr % 10000) == 0:
                    results_store.commit()
                    logger.info("Committing, %d records made", cntr)

            worker_pool.terminate()
        else:
            logger.debug("Building theoretical sequences sequentially")
            for row in self.ms1_results_reader:
                res = task_fn(row)
                for item in res:
                    results_store[cntr] = item
                    cntr += 1
                    if (cntr % 10000) == 0:
                        results_store.commit()
                        logger.info("Committing, %d records made", cntr)
        results_store.commit()
        results_store.close()
        return results_store.filename
