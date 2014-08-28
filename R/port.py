import pandas as pd
import json


json_serialized_columns = ["Oxonium_ions", "Stub_ions",
                           "b_ion_coverage", "b_ions_with_HexNAc",
                           "y_ion_coverage", "y_ions_with_HexNAc"]


def prepare_model_file(path):
    data = pd.read_csv(path)
    if "numOxIons" not in data:
        data["numOxIons"] = map(
            lambda x: len(json.loads(x)), data["Oxonium_ions"])
    if "numStubs" not in data:
        data["numStubs"] = map(
            lambda x: len(json.loads(x)), data["Stub_ions"])
    if "percent_b_ion_coverage" not in data:

    if "percent_y_ion_coverage" not in data:

    if "percent_b_ion_with_HexNAc_coverage" not in data:

    if "percent_y_ion_with_HexNAc_coverage" not in data:
