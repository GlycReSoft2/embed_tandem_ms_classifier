import csv

from glycresoft_ms2_classification.prediction_tools import convert_csv_to_nested
from glycresoft_ms2_classification.utils import try_get_outfile


def main(data_path, metadata_source=None, sites_file=None):
    prediction_results = convert_csv_to_nested(data_path)
    if metadata_source is not None:
        reader = csv.DictReader(open(metadata_source))
        glycans = get_glycan_identities(reader.fieldnames)
        prediction_results.metadata["glycan_identities"] = glycans
    if sites_file is not None:
        site_list = [line.strip() for line in open(sites_file, "r")]
        site_list = list(map(int, site_list))
        prediction_results.metadata["site_list"] = site_list
    outfile = try_get_outfile(data_path, "json")
    prediction_results.serialize(outfile)
    print(prediction_results.metadata)
    return outfile


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

if __name__ == '__main__':
    import argparse
    app = argparse.ArgumentParser()
    app.add_argument("csv_file")
    app.add_argument("-m", "--metadata", default=None)
    app.add_argument("-s", "--site-file", type=str, help="A file listing each position along\
     the peptide where glycosylation is predicted")
    args = app.parse_args()
    res = main(args.csv_file, args.metadata)
    print(res)
