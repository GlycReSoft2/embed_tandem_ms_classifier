import itertools
import json
import re
import sqlite3
import logging
from copy import deepcopy

from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


from ..structure.composition import composition_to_mass
from ..structure.sequence import Sequence
from ..structure.stub_glycopeptides import StubGlycopeptide

PROTON = composition_to_mass("p")
db_logger = logging.getLogger(__name__)


def neutral_mass(mz, z):
    return (mz * z) - (z * PROTON)


def mass_charge_ratio(neutral_mass, z):
    return (neutral_mass + (z * PROTON)) / z


class MSMSSqlDB(object):
    def __init__(self, connection_string=":memory:"):
        self.connection_string = connection_string
        self.connection = sqlite3.connect(connection_string)
        self.connection.row_factory = sqlite3.Row
        self.cursor = self.connection.cursor

    def init_schema(self):
        self.connection.executescript(ObservedPrecursorSpectrum.sql_schema())
        self.connection.executescript(ObservedTandemSpectrum.sql_schema())
        self.connection.commit()

    def load_data(self, precursors_list):
        for precursor in precursors_list:
            for stmt in precursor.to_sql():
                try:
                    self.connection.execute(stmt)
                except:
                    print(stmt)
                    raise
        self.connection.commit()

    def __getitem__(self, scan_id):
        results = []
        for row in self.execute('''select neutral_mass, charge,
             ObservedPrecursorSpectrum.precursor_id, other_data from ObservedPrecursorSpectrum
             join Scans on ObservedPrecursorSpectrum.precursor_id =
             Scans.precursor_id where scan_id = {0};'''.format(scan_id)):
            results.append(ObservedPrecursorSpectrum.from_sql(row, self))
        return results

    def __iter__(self):
        for row in self.execute("select * from ObservedPrecursorSpectrum;"):
            yield ObservedPrecursorSpectrum.from_sql(row, self)

    def execute(self, *args, **kwargs):
        return self.connection.execute(*args, **kwargs)

    def executemany(self, *args, **kwargs):
        return self.connection.executemany(*args, **kwargs)

    def executescript(self, *args, **kwargs):
        return self.connection.executescript(*args, **kwargs)

    def commit(self):
        self.connection.commit()

    def rollback(self):
        self.connection.rollback()

    def ppm_match_tolerance_search(self, mass, tolerance, target_table="ObservedPrecursorSpectrum",
                                   precursor_id=None, mass_shift=0):
        if mass_shift != 0:
            mass_expr = "(neutral_mass + {0})".format(mass_shift)
        else:
            mass_expr = "neutral_mass"
        query_str = "select * from {target_table} where\
         abs(({mass_expr} - {search_mass})/({mass_expr})) < {tolerance}"
        if precursor_id is not None:
            query_str += " and precursor_id={precursor_id}"
        results = self.execute(query_str.format(target_table=target_table, mass_expr=mass_expr,
                                                search_mass=mass,
                                                tolerance=tolerance, precursor_id=precursor_id))
        for result in results:
            yield result


class Scan(object):
    def __init__(self, id, z, mz):
        self.scan_id = id
        self.z = z
        self.mz = mz

    def __repr__(self):
        rep = "Scan({scan_id} {mz}|{z})".format(**self.__dict__)
        return rep


class ObservedPrecursorSpectrum(object):
    __table_name__ = "ObservedPrecursorSpectrum"

    @classmethod
    def sql_schema(cls):
        return '''
drop table if exists  ObservedPrecursorSpectrum;
create table ObservedPrecursorSpectrum (
    precursor_id integer unique primary key not null,
    charge integer,
    neutral_mass float,
    other_data text);

drop table if exists Scans;
create table Scans(
    scan_id integer unique primary key not null,
    mz float,
    z integer,
    precursor_id integer,
    foreign key(precursor_id) references ObservedPrecursorSpectrum(precursor_id));
'''

    def __init__(self, scan_data, scan_ids, z, neutral_mass, tandem_data, **data):
        self.charge = z
        self.scans = scan_data
        self.scan_ids = scan_ids
        self.neutral_mass = neutral_mass
        self.tandem_data = tandem_data
        self.other_data = data
        self._iterkey = None

    @property
    def id(self):
        return self._iterkey

    @id.setter
    def id(self, value):
        self._iterkey = value

    def __repr__(self):
        return "<Scans: ({scans}), Neutral Mass: {neutral_mass}, Tandem Spectra: {num_tandem}>".format(
            num_tandem=len(self.tandem_data),
            **self.__dict__)

    def __ge__(self, o):
        return self.neutral_mass >= o

    def __gt__(self, o):
        return self.neutral_mass > o

    def __lt__(self, o):
        return self.neutral_mass < o

    def __le__(self, o):
        return self.neutral_mass <= o

    def to_json(self):
        json_dict = deepcopy(self.__dict__)
        json_dict['tandem_data'] = [t.to_json()
                                    for t in json_dict['tandem_data']]
        return json_dict

    def to_sql(self):
        stmt = '''insert into {table} (precursor_id, neutral_mass, charge, other_data)
        VALUES ({id}, {neutral_mass}, {charge}, '{other_data}');'''.format(
            table="ObservedPrecursorSpectrum",
            id=self._iterkey, neutral_mass=self.neutral_mass, charge=self.charge,
            other_data=json.dumps(self.other_data))
        yield stmt

        for tandem in self.tandem_data:
            yield (tandem.to_sql(self._iterkey))

        for scan in self.scans:
            insert_stmt = '''insert into {table} (scan_id, mz, z, precursor_id)
            values ({id}, {mz}, {z},{precursor_id});'''.format(
                table="Scans", precursor_id=self._iterkey,  **scan)
            yield (insert_stmt)

    @classmethod
    def from_sql(cls, row, cursor):
        neutral_mass = row['neutral_mass']
        charge = row['charge']
        _iterkey = row['precursor_id']
        other_data = json.loads(row['other_data'])
        scans_rows = cursor.execute("select * from Scans where Scans.precursor_id={0};".format(_iterkey))
        scans = map(dict, scans_rows)
        scan_ids = [scan.pop('scan_id') for scan in scans]
        for i, scan in enumerate(scans):
            scan['id'] = scan_ids[i]
            scan.pop("precursor_id")
        tandem_rows = cursor.execute("select * from ObservedTandemSpectrum where\
         ObservedTandemSpectrum.precursor_id={0};".format(_iterkey))
        tandem_spectra = [ObservedTandemSpectrum.from_sql(row, cursor) for row in tandem_rows]
        instance = cls(scans, scan_ids, charge, neutral_mass, tandem_spectra, **other_data)
        instance._iterkey = _iterkey
        return instance


class ObservedTandemSpectrum(object):
    __table_name__ = "ObservedTandemSpectrum"

    @classmethod
    def sql_schema(cls):
        return '''
drop table if exists ObservedTandemSpectrum;
create table ObservedTandemSpectrum (
    tandem_id integer unique primary key not null,
    intensity float,
    charge integer,
    neutral_mass float,
    annotation text,
    other_data text,
    precursor_id integer,
    foreign key(precursor_id) references ObservedPrecursorSpectrum(precursor_id));
'''

    def __init__(self, mass, z, intensity, id=None, annotation=None, **data):
        if annotation is None:
            annotation = data.pop("annotation", [])
        self.charge = z
        self.mass = mass  # backwards compatibility
        self.neutral_mass = mass
        self.intensity = intensity
        self.id = data.pop("id", None) or id
        self.annotation = annotation
        self.other_data = data

    def __repr__(self):
        return "<ObservedTandemSpectra {neutral_mass}, {charge}, {intensity}>".format(**self.__dict__)

    def __ge__(self, o):
        return self.neutral_mass >= o

    def __gt__(self, o):
        return self.neutral_mass > o

    def __lt__(self, o):
        return self.neutral_mass < o

    def __le__(self, o):
        return self.neutral_mass <= o

    def __rge__(self, o):
        return self.neutral_mass <= o

    def __rgt__(self, o):
        return self.neutral_mass < o

    def __rlt__(self, o):
        return self.neutral_mass > o

    def __rle__(self, o):
        return self.neutral_mass >= o

    def to_json(self):
        return self.__dict__

    def to_sql(self, precursor_id=None):
        return '''insert into {table} (tandem_id, intensity, charge, neutral_mass, annotation, other_data, precursor_id)
        VALUES ({id}, {intensity}, {charge}, {neutral_mass}, '{annotation}', '{other_data}', {precursor_id});'''.format(
            table="ObservedTandemSpectrum",
            id=self.id, charge=self.charge, neutral_mass=self.neutral_mass, intensity=self.intensity,
            other_data=json.dumps(self.other_data), precursor_id=precursor_id, annotation=json.dumps(self.annotation)
        )

    @classmethod
    def from_sql(cls, row, cursor):
        mass = row['neutral_mass']
        charge = row['charge']
        other_data = json.loads(row['other_data'])
        annotation = json.loads(row["annotation"])
        intensity = row['intensity']
        id = row['tandem_id']
        instance = cls(mass, charge, intensity, id=id, annotation=annotation, **other_data)
        return instance


def extract_annotations(*spectra):
    return [dict(s.annotation) for s in spectra]


def matched_spectra(scan_ids, neutral_mass, ppm_error, match_key, **data):
    return dict(scan_ids=scan_ids,
                neutral_mass=neutral_mass,
                ppm_error=ppm_error,
                key=match_key,
                other_data=data)


def label_formatter(label):
    label = str(label).strip()
    if len(label) > 12:
        label = '\n'.join(
            s for s in re.split(r'([^\+-]+[\+-])', label) if len(s) != 0)

    return label, label.count("\n")


def plot_observed_spectra(spectra, title=None, label_formatter_fn=label_formatter,
                          include_target_ion_ladder=True, label_target_ion_ladder=False,
                          annotation_source=None):
    mzs = np.array(
        [mass_charge_ratio(spec.neutral_mass, spec.charge) for spec in spectra])
    intensities = np.array([spec.intensity for spec in spectra])
    charge = np.array([spec.charge for spec in spectra])
    fig = plt.figure()
    plt.xlabel("m/z")
    plt.ylabel("Relative Intensity")
    plt.bar(mzs, intensities, width=0.1, edgecolor='black')

    uppermost_text = 0
    label_map = extract_annotations(*spectra) if annotation_source is None else None
    if label_map is not None:
        track_colors = {}
        color_iter = itertools.cycle(colors)
        label_offset = itertools.cycle(label_shift)

        def get_color(key):
            try:
                return track_colors[key]
            except:
                track_colors[key] = color_iter.next()
                return track_colors[key]

        all_keys = {k for s in label_map for k in s.keys()}
        # Ensure that each key is associated with a color in case we have only
        # common ions
        map(get_color, all_keys)

        arrowprops = dict(width=0.01, headwidth=0.01, alpha=0.2)
        need_text_offset_reset = True
        for i, annots in enumerate(label_map):
            xy = (mzs[i], intensities[i])
            chg = charge[i]
            if need_text_offset_reset:
                offset_counter = label_offset.next()
                need_text_offset_reset = False

            all_same = len(set(map(str, annots.values()))) == 1 and len(
                annots) == len(all_keys) and (len(all_keys) != 1)
            if all_same:
                fmt_label, lines = label_formatter_fn(annots.values()[0])
                if chg > 1:
                    fmt_label = r"{1}$\tt{{^{{+{0}}}}}$".format(chg, fmt_label)
                plt.annotate(fmt_label, xy, arrowprops=arrowprops,
                             xytext=(-15, 10 + offset_counter + (lines * 12)),
                             textcoords='offset points',
                             color='black', fontweight=600,
                             va="top", ha="left")
                need_text_offset_reset = True
                offset_counter += (lines + 1) * 14
                uppermost_text = max(uppermost_text, offset_counter)
            else:
                for track, text in annots.items():
                    fmt_label, lines = label_formatter_fn(text)
                    if chg > 1:
                        fmt_label = r"{1}$\tt{{^{{+{0}}}}}$".format(chg, fmt_label)
                    plt.annotate(fmt_label, xy, arrowprops=arrowprops,
                                 xytext=(-15, 10 + offset_counter +
                                         (lines * 12)),
                                 textcoords='offset points',
                                 color=get_color(track),
                                 va="bottom", ha="left")
                    need_text_offset_reset = True
                    offset_counter += (lines + 1) * 14
                    uppermost_text = max(uppermost_text, offset_counter)
        # Generate the reference peaks expected from the matched precursors to
        # provide a visual guide for evaluating a match.
        if include_target_ion_ladder:
            max_intensity = max(intensities)
            for target_seq in all_keys:
                # Parse the key which is the sequence into an object and generate the fragments
                # as a flat list.
                seq_obj = Sequence(target_seq)
                b_frags = list(
                    itertools.chain.from_iterable(seq_obj.get_fragments('b')))
                y_frags = list(
                    itertools.chain.from_iterable(seq_obj.get_fragments('y')))

                stub_glycopeptides = StubGlycopeptide.from_sequence(seq_obj).get_stubs()

                # No charge information in fragment, use 1?
                b_frag_mz = [mass_charge_ratio(f.mass, 1) for f in b_frags]
                y_frag_mz = [mass_charge_ratio(f.mass, 1) for f in y_frags]
                stub_mz = [mass_charge_ratio(f["mass"], 1) for f in stub_glycopeptides]

                # Draw the bars for these reference ions. Make them visible but
                # don't wash out the actual observed ions
                plt.bar(b_frag_mz, [max_intensity] * len(b_frag_mz), width=0.01,
                        edgecolor=get_color(target_seq), alpha=0.1)
                plt.bar(y_frag_mz, [max_intensity] * len(y_frag_mz), width=0.01,
                        edgecolor=get_color(target_seq), alpha=0.1)
                plt.bar(stub_mz, [max_intensity] * len(stub_mz), width=0.01,
                        edgecolor=get_color(target_seq), alpha=0.1)
                if label_target_ion_ladder:
                    layer = 1
                    layer_size = 20
                    label_offset = itertools.cycle(label_shift_small)
                    for i, frag in enumerate(b_frags):
                        offset_counter = label_offset.next() + \
                            (layer * layer_size)
                        xy = (b_frag_mz[i], max_intensity)
                        label = frag.get_fragment_name()
                        fmt_label, lines = label_formatter_fn(label)
                        if len(fmt_label) > 4:
                            lines += 1
                        plt.annotate(fmt_label, xy, arrowprops=arrowprops,
                                     xytext=(
                                         0, 0 + offset_counter + (lines * 12)),
                                     textcoords='offset points',
                                     color=get_color(target_seq),
                                     va="bottom", ha="left")
                        uppermost_text = max(
                            uppermost_text, offset_counter + (lines * 12))

                    for i, frag in enumerate(y_frags):
                        offset_counter = label_offset.next() + \
                            (layer * layer_size)
                        xy = (y_frag_mz[i], max_intensity)
                        label = frag.get_fragment_name()
                        fmt_label, lines = label_formatter_fn(label)

                        plt.annotate(fmt_label, xy, arrowprops=arrowprops,
                                     xytext=(0, 0 + offset_counter +
                                             (lines * 12)),
                                     textcoords='offset points',
                                     color=get_color(target_seq),
                                     va="bottom", ha="left")
                        uppermost_text = max(uppermost_text, offset_counter)

                    for i, frag in enumerate(stub_glycopeptides):
                        offset_counter = label_offset.next() + \
                            (layer * layer_size)
                        xy = (stub_mz[i], max_intensity)
                        label = frag["key"]
                        fmt_label, lines = label_formatter_fn(label)

                        plt.annotate(fmt_label, xy, arrowprops=arrowprops,
                                     xytext=(
                                         0, 0 + offset_counter + (lines * 12)),
                                     textcoords='offset points',
                                     color=get_color(target_seq),
                                     va="bottom", ha="left")
                        uppermost_text = max(
                            uppermost_text, offset_counter + (lines * 12))
                    layer += 1

        plt.figlegend(map(lambda x: mpatches.Patch(color=x), track_colors.values()),
                      track_colors.keys(), "upper right", fontsize=12)
        plt.ylim(0, (max(intensities) + uppermost_text) * 1.25)
        plt.xlim(0, (max(mzs)) * 1.15)
    fig.set_size_inches(12, 10)
    return fig

label_shift = [60, 0, 40, 0, 20]
label_shift_small = [60, 0, 40, 0, 20]
colors = ['red', 'blue', 'orange', 'green', 'teal', 'darkcyan']
