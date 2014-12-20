from glycresoft_ms2_classification import error_code_interface
from glycresoft_ms2_classification.utils import csv_to_json
from StringIO import StringIO


def error_code_map():
    buffer_obj = StringIO()
    error_code_interface.ErrorCodingMeta.build_error_code_map(buffer_obj)
    print(buffer_obj.getvalue())


def convert_csv_to_json():
    csv_to_json()
