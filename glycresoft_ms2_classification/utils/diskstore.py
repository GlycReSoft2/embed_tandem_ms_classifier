import shelve
import os
import random
import atexit


def random_file_name(length=22):
    name = ''.join([chr(i) for i in random.sample(ascii_range, length)])
    return name
ascii_range = range(48, 58) + range(65, 91) + range(97, 123)

ALWAYS_NEW = 0


class DiskStore(object):

    def __init__(self, index_file_path=None, persistence_mode=ALWAYS_NEW):
        self.index_file_path = index_file_path
        self.persistence_mode = persistence_mode

        if persistence_mode == ALWAYS_NEW:
            if os.path.exists(index_file_path) and persistence_mode == ALWAYS_NEW:
                os.remove(index_file_path)
            atexit.register(lambda: os.remove(index_file_path))

        self.store = shelve.open(index_file_path)
        try:
            self.keys = self.store.iterkeys
            self.values = self.store.itervalues
            self.items = self.store.iteritems
        except:
            self.keys = self.store.keys
            self.values = self.store.values
            self.items = self.store.items

    def __getitem__(self, key):
        return self.store[key]

    def __setitem__(self, key, value):
        self.store[key] = value

    def __iter__(self):
        return iter(self.store)

    def __contains__(self, key):
        return key in self.store

    def popitem(self):
        return self.store.popitem()

    def close(self):
        self.store.close()

    def update(self, **kwargs):
        self.store.update(**kwargs)

    def refresh(self):
        self.close()
        os.remove(self.index_file_path)
        self.store = shelve.open(self.index_file_path)

    #Unused. Does not work with immutables
    def mixin_writeback_handle(self, key, obj):
        setattr(obj, "__backer", self)
        setattr(obj, "__backing_key", key)

        def writeback():
            obj.__backer[obj.__backing_key] = obj
        obj.writeback = writeback
