from ..structure.composition import composition_to_mass

PROTON = composition_to_mass("p")


def neutral_mass(mz, z):
    return (mz * z) - (z * PROTON)


class ObservedPrecursorSpectra(object):

    def __init__(self, scan_data, scan_ids, neutral_mass, tandem_data, **data):
        self.scans = scan_data
        self.scan_ids = scan_ids
        self.neutral_mass = neutral_mass
        self.tandem_data = tandem_data
        self.other_data = data
        self._iterkey = None

    def __repr__(self):
        return "<Scans: ({scans[0]}), Neutral Mass: {neutral_mass}, Tandem Spectra: {num_tandem}>".format(
            num_tandem=len(self.tandem_data),
            min_scan=self.scans[-1],
            **self.__dict__)

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


class ObservedTandemSpectra(object):

    def __init__(self, mass, z, intensity, **data):
        self.charge = z
        self.mass = mass
        self.intensity = intensity
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


def matched_spectra(scan_ids, neutral_mass, ppm_error, match_key, **data):
    return dict(scan_ids=scan_ids,
                neutral_mass=neutral_mass,
                ppm_error=ppm_error,
                key=match_key,
                other_data=data)
