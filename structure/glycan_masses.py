import csv
import os


mammalian_glycomedb_nlinked_path = os.path.join(os.path.dirname(__file__), "data", "Mammalian_GlycomeDB_NLinked.csv")


def load_from_file(path_to_file=mammalian_glycomedb_nlinked_path):
    with open(path_to_file) as fh:
        glycan_reader = csv.DictReader(fh)
        glycan_masses = list(set([("%0.3f" % float(g["Molecular Weight"]), g["Compositions"]) for g in glycan_reader]))
        return glycan_masses
