import multiprocessing
import logging
from Queue import Empty as QueueEmptyError
import sqlitedict


HALT = 1
logger = logging.getLogger("MSMS Merger")


def worker_task(annotation_store_path, match_ions_store_path, counter_store_path, data_queue, comm_line):
    annotation_store = sqlitedict.open(annotation_store_path, journal_mode="OFF")
    match_ions_store = sqlitedict.open(match_ions_store_path, journal_mode="OFF")
    did_match_counter = sqlitedict.open(counter_store_path, journal_mode="OFF")
    cntr = 0
    running = True
    while running:
        has_signal = comm_line.poll(0.01)
        if has_signal:
            signal = comm_line.recv()
            logger.debug("Signal found: %r", signal)
            if signal == HALT:
                running = False
        logger.debug("Processing Queue %d", cntr)
        while not data_queue.empty():
            try:
                matches, counter, annotater = data_queue.get(timeout=30)
                for m in matches:
                    match_ions_store[cntr] = m
                    cntr += 1
                increment_counter(did_match_counter, counter)
                combine_annotations(annotation_store, annotater)
            except QueueEmptyError:
                break
    annotation_store.commit()
    match_ions_store.commit()
    did_match_counter.commit()

    return did_match_counter


def combine_annotations(accumulator, collected):
    for k, v in collected.items():
        current = accumulator.get(k, [])
        current.extend(v)
        accumulator[k] = current


def increment_counter(accumulator, counter):
    for k, v in counter.items():
        current = accumulator.get(k, 0)
        accumulator[k] = current + v


class Merger(object):
    def __init__(self, db_path):
        self.annotation_store = ("./{}_.annotation_store".format(db_path))
        self.match_ions_store = ("./{}_.match_ions_store".format(db_path))
        self.counter_store = ("./{}_.counter_store".format(db_path))
        self.queue_manager = multiprocessing.Manager()
        self.data_queue = self.queue_manager.Queue()
        self.worker_process = None
        self.to_worker_comm_pipe, self.from_worker_comm_pipe = multiprocessing.Pipe()

    def start_worker(self):
        if self.worker_process is None or not self.worker_process.is_alive():
            arg_tuple = (self.annotation_store, self.match_ions_store, self.counter_store,
                         self.data_queue, self.from_worker_comm_pipe)
            self.worker_process = multiprocessing.Process(
                target=worker_task, args=arg_tuple, name="msms_merger")
            self.worker_process.start()

    def stop_worker(self):
        if self.worker_process.is_alive():
            logger.debug("Sending HALT signal")
            try:
                self.to_worker_comm_pipe.send(HALT)
                del self.worker_process
            except EOFError, e:
                print e

    def open_annotation_store(self):
        return sqlitedict.open(self.annotation_store)

    def open_match_ions_store(self):
        return sqlitedict.open(self.match_ions_store)

    def open_counter_store(self):
        return sqlitedict.open(self.counter_store)

    def sync_matches_db(self, other_db_file, table_name="matched_ions"):
        other_db = sqlitedict.open(other_db_file, table_name, journal_mode="OFF")
        matches = sqlitedict.open(self.match_ions_store)
        for k, v in matches.iteritems():
            other_db[k] = v
