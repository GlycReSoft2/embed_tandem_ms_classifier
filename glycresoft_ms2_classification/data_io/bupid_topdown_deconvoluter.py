import yaml
import itertools

from . import ObservedPrecursorSpectra
from . import ObservedTandemSpectra


MAX_CHARGE_STATE = 4


class BUPIDYamlParser(object):

    def __init__(self, file_path):
        self.file_path = file_path
        stream = open(file_path, 'r')
        try:
            loader = yaml.CLoader(stream)
        except:
            loader = yaml.Loader(stream)
        raw_data = (loader.get_data())
        self.data = dict()
        self._build_spectra(raw_data)

    def __iter__(self):
        return iter(self.data.items())

    def _build_spectra(self, raw_data):
        for tandem_ms_ind, peak_data in enumerate(raw_data['peaks']):
            scan_id_range = [scan["id"] for scan in peak_data["scans"]]
            # Treat the first scan as representative
            precursor = peak_data["scans"][0]
            precursor_mz = precursor["mz"]
            precursor_charge = precursor["z"]
            precursor_neutral_mass = precursor_mz / float(precursor_charge)
            tandem_data = [ObservedTandemSpectra(*ion) for
                           ion in itertools.izip(
                               peak_data["mass"], peak_data["z"], peak_data["intensity"])
                           if ion[1] <= MAX_CHARGE_STATE]

            observed_spectra = ObservedPrecursorSpectra(peak_data["scans"],
                                                        scan_id_range,
                                                        precursor_neutral_mass,
                                                        tandem_data)
            observed_spectra._iterkey = tandem_ms_ind
            self.data[scan_id_range[0]] = observed_spectra