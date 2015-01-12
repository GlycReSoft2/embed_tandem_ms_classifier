import pandas as pd
import json

from types import MethodType
from functools import partial


class MDFrameJSONEncoder(json.JSONEncoder):
    ''' '''
    def default(self, o):
        if isinstance(o, pd.DataFrame):
            return json.loads(o.to_json(orient="records"))
        else:
            return json.JSONEncoder.default(self, o)


class IndirectIndexer(object):
    '''A stubby front-end to mimic the DataFrame indexer components like .ix[],
    .loc[], and such. Since these use __getitem__ and __setitem__ to do their
    work rather than direct function calls, must use an object'''
    def __init__(self, wrapped_indexer_getter, wrapped_indexer_setter):
        self.wrapped_indexer_getter = wrapped_indexer_getter
        self.wrapped_indexer_setter = wrapped_indexer_setter

    def __getitem__(self, inds):
        return self.wrapped_indexer_getter(inds)

    def __setitem__(self, inds, value):
        self.wrapped_indexer_setter(inds, value)


class MDFrame(object):

    def __init__(self, metadata, data, **kwargs):
        self.metadata = metadata
        self.data = data

    def schema_upgrade(self, version):
        pass

    def __getitem__(self, inds):
        return self.predictions[inds]

    def __setitem__(self, inds, value):
        self.predictions[inds] = value

    def __getattr__(self, name):
        try:
            return object.__getattribute__(self, name)
        except:
            try:
                out = getattr(object.__getattribute__(self, "data"), name)
                if isinstance(out, MethodType):
                    out = partial(self._framewrapper_fn, name)
                return out
            except:
                raise AttributeError(
                    "{clstype} has no attribute {attr}".format(clstype=type(self).__name__, attr=name))

    # def __setattr__(self, name, value):
    #     try:
    #         object.__getattribute__(self, name)
    #         object.__setattr__(self, name, value)
    #     except:
    #         self[name] = value

    @classmethod
    def deserialize(cls, data_buffer, data_key, metadata_key="metadata", *args, **kwargs):
        if isinstance(data_buffer, basestring):
            data_buffer = open(data_buffer, "rb")
        loose_data = json.load(data_buffer)
        metadata = loose_data[metadata_key]
        try:
            content = pd.DataFrame(loose_data[data_key])
        except KeyError:
            print("Data Key %s not found. Keys: %s" % (data_key, loose_data.keys()))
            raise
        instance = cls(metadata, content, **kwargs)
        return instance

    def _framewrapper_fn(self, attr, *args, **kwargs):
        out = getattr(self.data, attr)(*args, **kwargs)
        if isinstance(out, pd.DataFrame):
            out = self.__class__(self.metadata, out)
        return out

    def __repr__(self):
        rep = "<MDFrame>\n{frame}".format(frame=repr(self.data))
        return rep

    def __iter__(self):
        return iter(self.data)

    def __dir__(self):
        return ["metadata"] + [attr for attr in dir(self.data) if "__" not in attr]

    def __contains__(self, name):
        return name in self.data

    def wrap_ix_get(self, index):
        df = self.data.ix[index]
        return self.__class__(self.metadata, df)

    def wrap_ix_set(self, index, value):
        self.data.ix[index] = value

    def wrap_loc_get(self, index):
        df = self.data.loc[index]
        return self.__class__(self.metadata, df)

    def wrap_loc_set(self, index, value):
        self.data.loc[index] = value

    def wrap_iloc_get(self, index):
        df = self.predictions.iloc[index]
        return self.__class__(self.metadata, df)

    def wrap_iloc_set(self, index, value):
        self.data.iloc[index] = value
