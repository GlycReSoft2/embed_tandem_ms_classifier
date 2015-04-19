import pkg_resources

import sqlitedict
sqlitedict.logger.setLevel("WARNING")

__all__ = ["entry_point", "structure", "prediction_tools", "proteomics",
           "match_ions", "theoretical_glycopeptide", "postprocess",
           "data_io", "utils"]

pkg_resources.declare_namespace("glycresoft_ms2_classification")
