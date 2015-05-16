import csv
import os
import re
import datetime
import multiprocessing
import logging

import functools
import itertools
import cPickle

from sqlitedict import SqliteDict

from glycresoft_ms2_classification.structure.modification import RestrictedModificationTable
from glycresoft_ms2_classification.structure.modification import ModificationTable
from glycresoft_ms2_classification.structure.sequence import Sequence
from glycresoft_ms2_classification.structure.sequence_space import SequenceSpace
from glycresoft_ms2_classification.structure.stub_glycopeptides import StubGlycopeptide
from glycresoft_ms2_classification.structure import constants

from glycresoft_ms2_classification.proteomics import get_enzyme, msdigest_xml_parser

logger = logging.getLogger("search_space_builder")
mod_pattern = re.compile(r'(\d+)(\w+)')
g_colon_prefix = "G:"


def parse_site_file(site_lists):
    site_list_map = {}
    last_name = None
    for line in site_lists:
        if len(line) == 0:
            continue
        elif ">" == line[0]:
            last_name = line[1:].strip()
        else:
            site_list_map[last_name] = map(int, re.findall(r"(\d+)\s+", line))
    return site_list_map


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
        num_missed_cleavages = int(kwargs.get("PeptideMissedCleavage#"))
        num_glycosylation_sites = int(kwargs.get("#ofGlycanAttachmentToPeptide"))
        mass_error = float(kwargs.get("PPM Error"))
        volume = float(kwargs.get("Total Volume"))
        start_pos = int(kwargs.get("StartAA"))
        end_pos = int(kwargs.get("EndAA"))
        protein_id = kwargs.get("ProteinID", None)
        glycan_composition_map = _get_glycan_counts(kwargs, glycan_identities)
        return cls(
            score=score, theoretical_precursor_mass=theoretical_precursor_mass,
            observed_precursor_mass=observed_precursor_mass, glycan_composition_str=glycan_composition_str,
            peptide_sequence=peptide_sequence, peptide_modifications=peptide_modifications,
            num_missed_cleavages=num_missed_cleavages, num_glycosylation_sites=num_glycosylation_sites,
            mass_error=mass_error, volume=volume, start_pos=start_pos, end_pos=end_pos,
            glycan_composition_map=glycan_composition_map, protein_id=protein_id)

    sql_schema = '''
    DROP TABLE IF EXISTS theoretical_search_space;
    CREATE TABLE theoretical_search_space(
        id INTEGER PRIMARY KEY,
        glycopeptide_identifier VARCHAR(200),
        structure TEXT
    );
    '''

    def __init__(self, score=None, theoretical_precursor_mass=None,
                 observed_precursor_mass=None, glycan_composition_str=None,
                 peptide_sequence=None, peptide_modifications=None, num_missed_cleavages=None,
                 num_glycosylation_sites=None, mass_error=None, volume=None, start_pos=None,
                 end_pos=None, glycan_composition_map=None, protein_id=None, fragments=None,
                 oxonium_ions=None, id=None, glycopeptide_identifier=None, modified_sequence=None):
        self.ms1_score = score
        self.theoretical_precursor_mass = theoretical_precursor_mass
        self.observed_precursor_mass = observed_precursor_mass
        self.glycan_composition_str = glycan_composition_str
        self.peptide_sequence = peptide_sequence
        self.peptide_modifications = peptide_modifications
        self.num_missed_cleavages = num_missed_cleavages
        self.num_glycosylation_sites = num_glycosylation_sites
        self.mass_error = mass_error
        self.volume = volume
        self.start_position = start_pos
        self.end_position = end_pos
        self.glycan_composition_map = glycan_composition_map
        self.protein_id = protein_id
        self.fragments = fragments
        self.oxonium_ions = oxonium_ions
        self.id = id
        self.glycopeptide_identifier = glycopeptide_identifier
        self.modified_sequence = modified_sequence

    def to_sql(self, id=None):
        stmt = '''INSERT OR  REPLACE INTO theoretical_search_space
                            (id, glycopeptide_identifier, structure)
                        VALUES
                            (?, ?, ?, ?, ?)'''
        self.id = id if id is not None else self.id
        return (stmt, (self.id, cPickle.dumps(self)))

    @classmethod
    def from_sql(cls, row, conn=None):
        return cPickle.loads(str(row['structure']))

    def as_dict(self):
        dict_form = {
            "MS1_Score": self.ms1_score,
            "Obs_Mass": self.observed_precursor_mass,
            "Calc_mass": self.theoretical_precursor_mass, "ppm_error": self.mass_error,
            "Peptide": self.peptide_sequence, "Peptide_mod": self.peptide_modifications,
            "Glycan": self.glycan_composition_str, "vol": self.volume, "glyco_sites": self.num_glycosylation_sites,
            "startAA": self.start_position, "endAA": self.end_position, "protein_id": self.protein_id
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
    b_ions_hexnac = []
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
                b_ions_hexnac.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})
            else:
                b_ions.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})

    y_type = fragments[1]  # seq.get_fragments('Y')
    y_ions = []
    y_ions_hexnac = []
    for y in y_type:
        for fm in y:
            key = fm.get_fragment_name()
            mass = fm.get_mass()
            golden_pairs = fm.golden_pairs
            if "HexNAc" in key:
                y_ions_hexnac.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})
            else:
                y_ions.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})

    pep_stubs = StubGlycopeptide(
        ms1_result.peptide_sequence,
        ms1_result.peptide_modifications,
        ms1_result.num_glycosylation_sites,
        ms1_result.glycan_composition_str)

    stub_ions = pep_stubs.get_stubs()
    oxonium_ions = pep_stubs.get_oxonium_ions()

    # Eventually, return this object instead of a dictionary
    ms1_result.glycopeptide_identifier = seq_mod + ms1_result.glycan_composition_str,
    ms1_result.modified_sequence = seq_mod
    ms1_result.oxonium_ions = oxonium_ions
    ms1_result.stub_ions = stub_ions
    ms1_result.bare_b_ions = b_ions
    ms1_result.bare_y_ions = y_ions
    ms1_result.b_ions_hexnac = b_ions_hexnac
    ms1_result.y_ions_hexnac = y_ions_hexnac

    fragments_dict = {
        "Seq_with_mod": seq_mod,
        "Glycopeptide_identifier": seq_mod + ms1_result.glycan_composition_str,
        "Oxonium_ions": oxonium_ions,
        "pep_stub_ions": stub_ions,
        "bare_b_ions": b_ions,
        "bare_y_ions": y_ions,
        "b_ions_with_HexNAc": b_ions_hexnac,
        "y_ions_with_HexNAc": y_ions_hexnac
    }

    fragments_dict.update(ms1_result.as_dict())

    return fragments_dict


def from_sequence(row, glycan_identity):
    """Convert an MS1GlycopeptideResult directly into its fragments
    with exact positions pre-specified on its :attr:`peptide_sequence`

    Parameters
    ----------
    row: dict
        A row from the input csv
    glycan_identity: list of str
        List of monosaccaride names

    Returns
    -------
    dict
    """
    ms1_result = MS1GlycopeptideResult.from_csvdict(glycan_identity, **row)
    if len(ms1_result.peptide_sequence) == 0:
        return []
    seq = Sequence(ms1_result.peptide_sequence)
    return [generate_fragments(seq, ms1_result)]


def process_predicted_ms1_ion(row, modification_table, site_list_map, glycan_identity):
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
    glycan_sites = set(site_list_map[ms1_result.protein_id]).intersection(
        range(ms1_result.start_position, ms1_result.end_position + 1))

    # No recorded sites, skip this component.
    if len(glycan_sites) == 0:
        return []

    # Adjust the glycan_sites to relative position
    glycan_sites = [x - ms1_result.start_position for x in glycan_sites]
    ss = get_search_space(
        ms1_result, glycan_sites, mod_list)
    seq_list = ss.get_theoretical_sequence(ms1_result.num_glycosylation_sites)
    fragments = [generate_fragments(seq, ms1_result)
                 for seq in seq_list]
    return fragments


class TheoreticalSearchSpace(object):
    '''
    Describe the process of generating all theoretical sequences and their fragments
    from an MS1 Results CSV, a collection of constant and variable peptide modifications,
    and a list of putative glycosylation sites.

    Uses an :class:`sqlitedict.SqliteDict` instance as a storage backend.

    Includes a single- and multi-process compatible implementation. The more processes used,
    the more memory must be allocated to buffer results.
    '''
    table_name = "theoretical_search_space"

    def __init__(self, ms1_results_file, db_file_name=None,
                 enzyme_info=None,
                 site_list=None,
                 constant_modifications=None,
                 variable_modifications=None, n_processes=4):
        if db_file_name is None:
            db_file_name = os.path.splitext(ms1_results_file)[0] + '.db'
        self.store = SqliteDict(db_file_name, tablename=self.table_name, journal_mode="OFF")
        self.metadata = SqliteDict(db_file_name, tablename='metadata', journal_mode="OFF")
        self.ms1_results_file = ms1_results_file

        try:
            site_list_map = parse_site_file(open(site_list))
        except IOError, e:
            if isinstance(site_list, dict):
                site_list_map = site_list
            else:
                raise e

        modification_table = RestrictedModificationTable.bootstrap(constant_modifications, variable_modifications)
        if constant_modifications is None and variable_modifications is None:
            modification_table = ModificationTable.bootstrap()
        self.modification_table = modification_table
        self.glycosylation_site_map = site_list_map
        self.ms1_results_reader = csv.DictReader(open(self.ms1_results_file))

        self.n_processes = n_processes

        self.glycan_identities = get_glycan_identities(self.ms1_results_reader.fieldnames)
        enzyme_info = map(get_enzyme, enzyme_info)
        tag = datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d-%H%M%S")
        self.metadata.update({
            "glycan_identities": self.glycan_identities,
            "enzyme_info": enzyme_info,
            "site_list_map": site_list_map,
            "constant_modification_list": constant_modifications,
            "variable_modification_list": variable_modifications,
            "ms1_output_file": ms1_results_file,
            "enzyme": enzyme_info,
            "tag": tag,
            "enable_partial_hexnac_match": constants.PARTIAL_HEXNAC_LOSS
        })
        self.metadata.commit()
        self.metadata.close()

    def prepare_task_fn(self):
        task_fn = functools.partial(process_predicted_ms1_ion, modification_table=self.modification_table,
                                    site_list_map=self.glycosylation_site_map, glycan_identity=self.glycan_identities)
        return task_fn

    def run(self):
        '''
        Execute the algorithm on :attr:`n_processes` processes
        '''
        task_fn = self.prepare_task_fn()
        results_store = self.store
        cntr = 0
        if self.n_processes > 1:
            worker_pool = multiprocessing.Pool(self.n_processes)
            logger.debug("Building theoretical sequences concurrently")
            for res in (itertools.chain.from_iterable(
                    worker_pool.imap(task_fn, self.ms1_results_reader, chunksize=500))):
                res['id'] = cntr
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
                    item['id'] = cntr
                    results_store[cntr] = item
                    cntr += 1
                    if (cntr % 10000) == 0:
                        results_store.commit()
                        logger.info("Committing, %d records made", cntr)
        results_store.commit()
        results_store.close()
        return results_store.filename


class ExactSearchSpace(TheoreticalSearchSpace):
    def __init__(self, ms1_results_file, db_file_name=None,
                 enzyme_info=None,
                 site_list=None,
                 constant_modifications=None,
                 variable_modifications=None, n_processes=4):

        if db_file_name is None:
            db_file_name = os.path.splitext(ms1_results_file)[0] + '.db'
        self.store = SqliteDict(db_file_name, tablename=self.table_name, journal_mode="OFF")
        self.metadata = SqliteDict(db_file_name, tablename='metadata', journal_mode="OFF")
        self.ms1_results_file = ms1_results_file

        try:
            site_list_map = parse_site_file(open(site_list))
        except IOError, e:
            if isinstance(site_list, dict):
                site_list_map = site_list
            else:
                raise e
        except TypeError:
            site_list_map = None

        self.ms1_results_reader = csv.DictReader(open(self.ms1_results_file))

        self.n_processes = n_processes

        self.glycan_identities = get_glycan_identities(self.ms1_results_reader.fieldnames)
        enzyme_info = map(get_enzyme, enzyme_info)
        tag = datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d-%H%M%S")
        self.metadata.update({
            "glycan_identities": self.glycan_identities,
            "enzyme_info": enzyme_info,
            "site_list_map": site_list_map,
            "constant_modification_list": constant_modifications,
            "variable_modification_list": variable_modifications,
            "ms1_output_file": ms1_results_file,
            "enzyme": enzyme_info,
            "tag": tag,
            "enable_partial_hexnac_match": constants.PARTIAL_HEXNAC_LOSS
        })
        self.metadata.commit()
        self.metadata.close()

    def prepare_task_fn(self):
        task_fn = functools.partial(from_sequence, glycan_identity=self.glycan_identities)
        return task_fn


class MS1ResultsFile(object):
    def __init__(self, file_path):
        self.file_path = file_path
        self.file_handle = open(file_path, 'rb')
        self.reader = csv.DictReader(self.file_handle)
        self.glycan_identities = get_glycan_identities(self.reader.fieldnames)

    def __iter__(self):
        for row in self.reader:
            yield MS1GlycopeptideResult.from_csvdict(self.glycan_identities, **row)
