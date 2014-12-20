import pkg_resources

__all__ = ["entry_point", "structure", "prediction_tools", "proteomics",
           "match_ions", "theoretical_glycopeptide", "postprocess",
           "data_io", "utils"]

try:
    from prediction_tools import prepare_model_file
    from prediction_tools import save_model_file
except:
    pass

pkg_resources.declare_namespace("glycresoft_ms2_classification")
