import logging
import sqlitedict

logger = logging.getLogger("disk_queue")
logger.setLevel("DEBUG")
sqlitedict.logger.setLevel('WARN')


class SqliteDiskQueue(object):
    def __init__(self, path=None):
        self.store = sqlitedict.SqliteDict(path, journal_mode="OFF")
        self.index = 0

    def append(self, obj):
        self.store[self.index] = obj
        self.index += 1

    def __iter__(self):
        for obj in self.store.itervalues():
            yield obj

    def __len__(self):
        return self.index
