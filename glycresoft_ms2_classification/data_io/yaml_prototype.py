import itertools
import yaml
import pandas as pd


def series_to_yaml_safe(series):
    index = series.index.to_native_types()
    values = series.values.tolist()

    return {i: v for i, v in itertools.izip(index, values)}


def frame_to_yaml_safe(frame):
    return {col: series_to_yaml_safe(series)
            for col, series in frame.iteritems()}


class YAMLFrame(yaml.YAMLObject):
    yaml_tag = u"PredictionResultsFrame"

    def __init__(self, metadata, **data):
        self.metadata = metadata
        self.frame = pd.DataFrame(data)

    def __getstate__(self):
        return {
            "metadata": self.metadata,
            "frame": frame_to_yaml_safe(self.frame)
        }

    def __setstate__(self, state):
        print(state)
        self.metadata = state['metadata']
        self.frame = pd.DataFrame(state['frame'])

    def __repr__(self):
        return yaml.dump(self)
    # def construct_yaml_object(self, node, cls):
    #     data = cls.__new__(cls)
    #     yield data
    #     if hasattr(data, '__setstate__'):
    #         state = self.construct_mapping(node, deep=True)
    #         data.__setstate__(state)
    #     else:
    #         state = self.construct_mapping(node)
    #         data.__dict__.update(state)
