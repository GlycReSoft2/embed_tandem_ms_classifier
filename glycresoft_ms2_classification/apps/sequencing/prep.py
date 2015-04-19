from glycresoft_ms2_classification import prediction_tools
from glycresoft_ms2_classification.ms import default_loader

predictions = prediction_tools.prepare_model_file("./ResultOf20131219_006_isos.scored_fdr.json").query("MS2_Score >= 0.5")
observed = default_loader("./20131219_006.mzML.20150125-173417matched.pkl").index_by_scan_ids()

spec = observed[predictions.iloc[2].scan_id_range[0]]
drop_mass = predictions.iloc[2].glycanMass
